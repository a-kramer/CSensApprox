# PHI

Approximating using the trapezoidal rule, with `DT=t-t0`:
```
I0(tf;yi) = I0 = identity(ny,ny)
I1(tf;ti) = integrate from ti to tf: A(t) I0(t;ti) dt ≈ 0.5 DT { A(tf)*I0 + A(ti)*I0 }
I2(tf;ti) = integrate from ti to tf: A(t) I1(t;ti) dt ≈ 0.5 DT { A(tf)*I1(tf;ti) + A(ti)*I1(ti;ti) }
                                                                                         ^^^^^^^^
																						    0.0
I3(tf;ti) = integrate from ti to tf: A(t)*0.5*DT*A(t)*I1(t;ti) dt ≈ [0.5*DT]^2 { A(tf)*A(tf)*I1(tf;ti) }
```
with 
```
s=0.5*DT 
I1 = s { A(tf) + A(ti) }
```
we sum the terms of PHI up:
```
PHI = 1 + I1 + s*A(tf)*I1(tf;ti) + s^2 * A(tf)^2 * I1(tf;ti) + [...]
```
in reverse order
```
PHI = (((s*A(tf) + 1)*s*A(tf) + 1)*s*A(tf) + 1)*I1(tf;ti) + I0
```
 
 
# S
 
```
S = PHI(k+1,k) * (S + 0.5*PHI(k,k)\B(k) + 0.5*PHI(k+1,k)\B(k+1))
```

with `PHI(k,k)=identity` (`1`) we re-write:

```
S = PHI(k+1,k) * (S + 0.5*(B(k) + PHI(k+1,k)\B(k+1))
S = PHI(k+1,k) * (S + 0.5*(B(k) + PHI(k,k+1)*B))
```

```
PHI(k+1,k)\B(k+1) = PHI(k,k+1)*B(k+1)
```


