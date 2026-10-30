"""Independent exact-Fraction reproduction; no imports from experiments.

Bounds each positive derivative series at the interval endpoints directly,
rather than using the certificate's fixed-binary interval Horner arithmetic.
Analytic conversion to Hessian/tensor bounds is the same mathematical input.
This independently verifies contraction with a coarser bound <0.997328724.
The separate fine-grid verifier proves the sharper <0.990737 bound.
"""
from fractions import Fraction as F
from math import factorial
from itertools import product
from collections import defaultdict, Counter
from pathlib import Path
import json

def require(condition):
    if not condition:
        raise RuntimeError('Independent rational or geometry check failed')

def add(a,b): return (a[0]+b[0], a[1]+b[1])
def neg(a): return (-a[1],-a[0])
def sub(a,b): return add(a,neg(b))
def mul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)
def scale(a,b): return mul(a,(F(b),F(b)))
def div(a,b):
    require(not b[0] <= 0 <= b[1])
    return mul(a,(1/b[1],1/b[0]))
def power(a,n):
    out=(F(1),F(1))
    for _ in range(n): out=mul(out,a)
    return out
def absup(a): return max(abs(x) for x in a)

R=F(12,5); DEG=18; CELLS=1024
coeff=[[F(1,4**(m+j)*factorial(m)*factorial(m+j+1)) for m in range(DEG+2)] for j in range(4)]
ratios=[R**2/(4*(DEG+2)*(DEG+j+3)) for j in range(4)]
require(all(0<=q<1 for q in ratios))
tails=[coeff[j][-1]*R**(2*(DEG+1))/(1-ratios[j]) for j in range(4)]
def zbound(j,s):
    # Positive coefficients imply exact endpoint extrema for the polynomial.
    lo=sum((coeff[j][m]*s[0]**m for m in range(DEG+1)),F(0))
    hi=sum((coeff[j][m]*s[1]**m for m in range(DEG+1)),F(0))+tails[j]
    return lo,hi

lam=F(1); third=F(0)
for k in range(CELLS):
    r=(R*k/CELLS,R*(k+1)/CELLS); ss=(r[0]**2,r[1]**2)
    z0,z1,z2,z3=[zbound(j,ss) for j in range(4)]
    p=div(z1,z0)
    q=sub(div(z2,z0),power(p,2))
    u=add(sub(div(z3,z0),scale(div(mul(z1,z2),power(z0,2)),3)),scale(power(p,3),2))
    trans=scale(p,2); radial=add(trans,scale(mul(ss,q),4))
    lam=min(lam,trans[0],radial[0])
    a=scale(mul(r,q),12); gamma=scale(mul(power(r,3),u),8)
    third=max(third,absup(add(a,gamma)),2*absup(a)/3)

# Deliberately coarser simple constants, checked against the entire cover.
L=F(137,1000); M=F(65,1000)
require(lam>=L and third<=M)
b=F(77,200); tau=F(1,4); delta=tau-L
mean_r=5*b; zz=(mean_r**2,mean_r**2)
mean=2*mean_r*zbound(1,zz)[1]/zbound(0,zz)[0]
w=3*b*b*delta
def exp_bound(x,n=40):
    term=total=F(1)
    for k in range(1,n+1):
        term=term*x/k; total+=term
    ratio=x/(n+2); require(0<=ratio<1)
    return total+term*x/(n+1)/(1-ratio)
ew=exp_bound(w)
require(F(289,1000)**2 > F(1,12))
s=tau+(ew-1)*ew*F(289,1000)
alpha=s*(6*b*mean+60*b*b*tau)+ew*(12*b*b*delta+60*b**3*M)/3
require(alpha<1)

require(alpha <= F(997328724, 10**9) < 1)
require(6*b<R)

# Independent unwrapped incidence enumeration, table data read without imports.
# AST literal data only; no source program execution.
import ast
def table(path,name):
    tree=ast.parse(Path(__file__).resolve().parent.parent.joinpath(path).read_text())
    return next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in n.targets))
colors=table('experiments/su2_color_geometry.py','COLOR_ROWS')
sigma=table('experiments/su2_forest_geometry.py','SIGMA_ROWS')
def color(e):
    x,d=e; return colors[tuple(v%2 for v in x)][d]
def shift(x,d): return tuple(v+(i==d) for i,v in enumerate(x))
def face(x,u,v): return ((x,u),(shift(x,u),v),(shift(x,v),u),(x,v))
stars=defaultdict(list); faces=0
for x in product(range(-2,5),repeat=4):
    for u in range(4):
        for v in range(u+1,4):
            edges=face(x,u,v)
            require(sorted(color(e) for e in edges)==[-1,0,1,2])
            eliminated=next(e for e in edges if color(e)==2)
            stars[eliminated].append(tuple(e for e in edges if color(e) in (0,1)))
            faces+=1
counts=Counter()
for x in product(range(2),repeat=4):
    require(colors[x][sigma[x][0]]==-1)
    for d in range(4):
        target=(x,d)
        if color(target) not in (0,1): continue
        occurrences=[]; same=cross=0
        for e,pairs in stars.items():
            if not any(target in pair for pair in pairs): continue
            require(len(pairs)==6)
            require(len({q for pair in pairs for q in pair})==12)
            own=next(pair for pair in pairs if target in pair)
            same+=1
            cross+=sum(len(pair) for pair in pairs if pair!=own)
            occurrences.extend(q for pair in pairs for q in pair if q!=target)
        require(same==6 and cross==60)
        require(len(occurrences)==66 and len(set(occurrences))==52)
        counts[(same,cross,len(set(occurrences)))]+=1
require(counts=={(6,60,52):32})

def lower_decimal(q,den=10**12): return str(F(q.numerator*den//q.denominator,den))
def upper_decimal(q,den=10**12): return str(F(-((-q.numerator*den)//q.denominator),den))
summary = {'independent_cells':CELLS,'independent_degree':DEG,
    'lambda_lower_outward':lower_decimal(lam),'third_upper_outward':upper_decimal(third),
    'simple_lambda_used':str(L),'simple_third_used':str(M),
    'independent_alpha_upper_outward':upper_decimal(alpha),
    'scope':'Independent coarser contraction check only; run verify.py for alpha<0.990737.',
    'certified_upper_threshold_rational':str(F(997328724,10**9)),
    'unwrapped_interior_links_checked':32,'unwrapped_faces_checked':faces,
    'unwrapped_counts':'each link: 6 same, 60 cross, 52 distinct neighbors',
    'floating_point_used':False}
rendered=json.dumps(summary,indent=2,sort_keys=True)+'\n'
Path(__file__).with_name('independent-summary.json').write_text(rendered,encoding='utf-8')
print(rendered,end='')
