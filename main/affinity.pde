class Affinity
{
  pt M1, M2, T; // Columns of a 2x3 matrix
  pt F; // fixed point
  
  // SVD
  pt U1, U2, VT1, VT2, S;
  float dc, z, l1r, l1i, l2r, l2i;
  
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
    float e11 = 0, e12 = 0, e21 = 0, e22 = 0;
    if (dc >= 0) {
      e11 = 0.5 * (exp(l1r * t) + exp(l2r * t)) + 0.5 * z / sqrt(dc) * (exp(l1r * t) - exp(l2r * t));
      e12 = 0.5 * M2.x / sqrt(dc) * (exp(l1r * t) - exp(l2r * t));
      e21 = 0.5 * M1.y / sqrt(dc) * (exp(l1r * t) - exp(l2r * t));
      e22 = 0.5 * (exp(l1r * t) + exp(l2r * t)) - 0.5 * z / sqrt(dc) * (exp(l1r * t) - exp(l2r * t));
    } else {
      e11 = exp(l1r * t) * cos(l1i * t) + z / sqrt(-dc) * exp(l1r * t) * sin(l1i * t);
      e12 = M2.x / sqrt(-dc) * exp(l1r * t) * sin(l1i * t);
      e21 = M1.y / sqrt(-dc) * exp(l1r * t) * sin(l1i * t);
      e22 = exp(l1r * t) * cos(l1i * t) - z / sqrt(-dc) * exp(l1r * t) * sin(l1i * t);
    }
    
    //vec FPnew = V(e11 * FP.x + e12 * FP.y, e21 * FP.x + e22 * FP.y);
    //if (debug) {
    //  println("t: " + t);
    //  println("E");
    //  println(e11 + ", " + e12);
    //  println(e21 + ", " + e22);
    //}
    pt Pnew = P(e11 * FP.x + e12 * FP.y + F.x, e21 * FP.x + e22 * FP.y + F.y);
    return Pnew;
  }
}
