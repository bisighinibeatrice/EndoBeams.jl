# EndoBeams: beam-to-beam contact 

Based on a point-to-point formulation. One new example: beam2beam.jl in /examples. 

Differences with respect to EndoBeams2.0:
- Remove use of sparse matrices 
- Use of LinearSolve.jl 

Still to be improved: 
- Beam-to-beam contact matrix computation 
- Canditates search 

----------------------------
## References
[1] C. Meier, W. A. Wall, A. Popp, A unified approach for beam-to-beam contact, Computer Methods in Applied Mechanics
and Engineering 315 (2017) 972-1010.
[2] P. Wriggers, G. Zavarise, On contact between three-dimensional beams undergoing large de
ections, Communications in
Numerical Methods in Engineering 13 (1997) 429-438.
[3] T. Otani, S.Wada, M. Tanaka, Modeling of endovascular coiling for cerebral aneurysms: Effects of friction on coil mechanical
behaviors, International Journal of Mechanical Sciences 166 (2020) 105206.
[4] C. Meier, A. Popp, W. A. Wall, Geometrically Exact Finite Element Formulations for Slender Beams: Kirchhoff-Love
Theory Versus Simo-Reissner Theory, 2017.
[5] M. Aguirre, S. Avril, An implicit 3D corotational formulation for frictional contact dynamics of beams against rigid surfaces
using discrete signed distance fields, Computer Methods in Applied Mechanics and Engineering 371 (2020) 113275.
[6] P. Wriggers, T. Lausen, Computational Contact Mechanics, 2008.