class Affinity
{
  pt M1, M2, T; // Columns of a 2x3 matrix
  pt F; // fixed point
  
  // SVD
  pt U1, U2, VT1, VT2, S;
  
  Affinity() {
    M1 = new pt(); M2 = new pt(); T = new pt(); F = new pt();
    U1 = new pt(); U2 = new pt(); S = new pt(); VT1 = new pt(); VT2 = new pt();
  }
  
  vec velocity(pt p) {
    return V(
      p.x * M1.x + p.y * M2.x + T.x,
      p.x * M1.y + p.y * M2.y + T.y);
  }
  
  pt apply(pt P, float t) {
    // Not correct
    vec FP = V(F, P);
    float sx = exp(t * S.x);
    float sy = exp(t * S.y);
    pt US1 = P(U1.x * sx, U1.y * sx);
    pt US2 = P(U2.x * sy, U2.y * sy);
    pt USVT1 = P(US1.x * VT1.x + US2.x * VT1.y, US1.y * VT1.x + US2.y * VT1.y);
    pt USVT2 = P(US1.x * VT2.x + US2.x * VT2.y, US1.y * VT2.x + US2.y * VT2.y);
    vec FPnew = V(USVT1.x * FP.x + USVT2.x * FP.y, USVT1.y * FP.x + USVT2.y * FP.y);
    return P(F, FPnew);
  }
}
