vec advectVector(vec V, pt P, float t, pt A, vec Av , pt B, vec Bv, pt C, vec Cv) {
   return V(
     advectPoint(P, t, A, Av, B, Bv, C, Cv),
     advectPoint(P(P, V), t, A, Av, B, Bv, C, Cv)
   ); 
}


pt advectPoint(pt P, float t, pt A, vec Av , pt B, vec Bv, pt C, vec Cv) {
  pt P1 = A, P2 = B, P3 = C;
  vec V1 = Av, V2 = Bv, V3 = Cv;
  pt M1 = new pt(), M2 = new pt(), T = new pt(), F = new pt();
  
  // Determinant
  float det = P1.x * P2.y + P2.x * P3.y + P3.x * P1.y - P3.x * P2.y - P2.x * P1.y - P1.x * P3.y;
  
  // Inverse Matrix
  float i11 = (P2.y - P3.y) / det;
  float i12 = (P3.x - P2.x) / det;
  float i13 = (P2.x * P3.y - P3.x * P2.y) / det;
  float i21 = (P3.y - P1.y) / det;
  float i22 = (P1.x - P3.x) / det;
  float i23 = (P3.x * P1.y - P3.y * P1.x) / det;
  float i31 = (P1.y - P2.y) / det;
  float i32 = (P2.x - P1.x) / det;
  float i33 = (P1.x * P2.y - P1.y * P2.x) / det;
  
  // Output
  M1.x = V1.x * i11 + V2.x * i21 + V3.x * i31;
  M2.x = V1.x * i12 + V2.x * i22 + V3.x * i32;
  T.x = V1.x * i13 + V2.x * i23 + V3.x * i33;
  M1.y = V1.y * i11 + V2.y * i21 + V3.y * i31;
  M2.y = V1.y * i12 + V2.y * i22 + V3.y * i32;
  T.y = V1.y * i13 + V2.y * i23 + V3.y * i33;
  
  det = M1.x * M2.y - M1.y * M2.x;
  i11 = M2.y / det;
  i12 = -M2.x / det;
  i21 = -M1.y / det;
  i22 = M1.x / det;
  F.x = -i11 * T.x - i12 * T.y;
  F.y = -i21 * T.x - i22 * T.y;
  
  // Matrix exp
  
  float dc, z;
  float l1r, l2r, l1i, l2i;
  
  dc = (M1.x + M2.y) * (M1.x + M2.y) / 4. - (M1.x * M2.y - M1.y * M2.x);
  z = (M1.x - M2.y) / 2.;
  if (dc >= 0) {
    l1r = (M1.x + M2.y) / 2. + sqrt(dc);
    l2r = (M1.x + M2.y) / 2. - sqrt(dc);
    l1i = 0;
    l2i = 0;
  } else {
    l1r = (M1.x + M2.y) / 2.;
    l2r = (M1.x + M2.y) / 2.;
    l1i = sqrt(-dc);
    l2i = -sqrt(-dc);
  }
  
  // Apply time t
  
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
    
    pt Pnew = P(e11 * FP.x + e12 * FP.y + F.x, e21 * FP.x + e22 * FP.y + F.y);
    return Pnew;
}

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

void compute_velocity_affinity(pt P1, pt P2, pt P3, vec V1, vec V2, vec V3, Affinity a) {
  // input positions: P1, P2, P3
  // input velocities: V1, V2, V3
  // output affine matrix: [M1, M2, T]
  // computes [V1 V2 V3] = [M1 M2 T] [P1 P2 P3//1 1 1]
  
  // Determinant
  float det = P1.x * P2.y + P2.x * P3.y + P3.x * P1.y - P3.x * P2.y - P2.x * P1.y - P1.x * P3.y;
  
  // Inverse Matrix
  float i11 = (P2.y - P3.y) / det;
  float i12 = (P3.x - P2.x) / det;
  float i13 = (P2.x * P3.y - P3.x * P2.y) / det;
  float i21 = (P3.y - P1.y) / det;
  float i22 = (P1.x - P3.x) / det;
  float i23 = (P3.x * P1.y - P3.y * P1.x) / det;
  float i31 = (P1.y - P2.y) / det;
  float i32 = (P2.x - P1.x) / det;
  float i33 = (P1.x * P2.y - P1.y * P2.x) / det;
  
  // Output
  a.M1.x = V1.x * i11 + V2.x * i21 + V3.x * i31;
  a.M2.x = V1.x * i12 + V2.x * i22 + V3.x * i32;
  a.T.x = V1.x * i13 + V2.x * i23 + V3.x * i33;
  a.M1.y = V1.y * i11 + V2.y * i21 + V3.y * i31;
  a.M2.y = V1.y * i12 + V2.y * i22 + V3.y * i32;
  a.T.y = V1.y * i13 + V2.y * i23 + V3.y * i33;
  
  if (debug) {
    println("Calculate Affinity");
    println("P1: (" + P1.x + ", " + P1.y + ")");
    println("P2: (" + P2.x + ", " + P2.y + ")");
    println("P3: (" + P3.x + ", " + P3.y + ")");
    println("V1: (" + V1.x + ", " + V1.y + ")");
    println("V2: (" + V2.x + ", " + V2.y + ")");
    println("V3: (" + V3.x + ", " + V3.y + ")");
    //println("Det: " + det);
    //println("Inv:");
    //println(i11 + " " + i12 + " " + i13);
    //println(i21 + " " + i22 + " " + i23);
    //println(i31 + " " + i32 + " " + i33);
    println("Affine:");
    println(a.M1.x + " " + a.M2.x + " " + a.T.x);
    println(a.M1.y + " " + a.M2.y + " " + a.T.y);
  }
  
  // Compute fixed point
  // F = (I-M)^(-1)T
  
  //det = (a.M1.x-1) * (a.M2.y-1) - a.M1.y * a.M2.x;
  //i11 = (a.M2.y-1) / det;
  //i12 = -a.M2.x / det;
  //i21 = -a.M1.y / det;
  //i22 = (a.M1.x - 1) / det;
  //a.F.x = -i11 * a.T.x - i12 * a.T.y;
  //a.F.y = -i21 * a.T.x - i22 * a.T.y;
  
  // Compute fixed point, modified version
  //F = (-M)^(-1)T
  
  det = a.M1.x * a.M2.y - a.M1.y * a.M2.x;
  i11 = a.M2.y / det;
  i12 = -a.M2.x / det;
  i21 = -a.M1.y / det;
  i22 = a.M1.x / det;
  a.F.x = -i11 * a.T.x - i12 * a.T.y;
  a.F.y = -i21 * a.T.x - i22 * a.T.y;
  
  //if (debug) {
  //  //println("Det: " + det);
  //  //println("Inv:");
  //  //println(i11 + " " + i12);
  //  //println(i21 + " " + i22);
  //  println("Fixed point: (" + a.F.x + ", " + a.F.y + ")");
  //}
  
  // Compute SVD
  // M = U * S * VT
  // https://lucidar.me/en/mathematics/singular-value-decomposition-of-a-2x2-matrix/
  
  //float Urx = 2 * a.M1.x * a.M1.y + 2 * a.M2.x * a.M2.y;
  //float Ury = a.M1.x * a.M1.x + a.M2.x * a.M2.x - a.M1.y * a.M1.y - a.M2.y * a.M2.y;
  //float Ur = 0.5 * atan2(Urx, Ury);
  //a.U1.x = cos(Ur);
  //a.U1.y = sin(Ur);
  //a.U2.x = -sin(Ur);
  //a.U2.y = cos(Ur);
  
  //float s1 = a.M1.x * a.M1.x + a.M2.x * a.M2.x + a.M1.y * a.M1.y + a.M2.y * a.M2.y;
  //float s2 = sqrt(Ury * Ury + Urx * Urx);
  //a.S.x = sqrt((s1 + s2) * 0.5);
  //a.S.y = sqrt((s1 - s2) * 0.5);
  
  //float Vrx = 2 * a.M1.x * a.M2.x + 2 * a.M1.y * a.M2.y;
  //float Vry = a.M1.x * a.M1.x - a.M2.x * a.M2.x + a.M1.y * a.M1.y - a.M2.y * a.M2.y;
  //float Vr = 0.5 * atan2(Vrx, Vry);
  //float s11 = (a.M1.x * cos(Ur) + a.M1.y * sin(Ur)) * cos(Vr) + (a.M2.x * cos(Ur) + a.M2.y * sin(Ur)) * sin(Vr);
  //float s22 = (a.M1.x * sin(Ur) - a.M1.y * cos(Ur)) * sin(Vr) + (-a.M2.x * sin(Ur) + a.M2.y * cos(Ur)) * cos(Vr);
  //float sgn11 = s11 / abs(s11);
  //float sgn22 = s22 / abs(s22);
  //a.VT1.x = cos(Vr) * sgn11;
  //a.VT1.y = sin(Vr) * sgn11;
  //a.VT2.x = -sin(Vr) * sgn22;
  //a.VT2.y = cos(Vr) * sgn22;
  
  // Matrix exponential in SAM paper
  a.dc = (a.M1.x + a.M2.y) * (a.M1.x + a.M2.y) / 4. - (a.M1.x * a.M2.y - a.M1.y * a.M2.x);
  a.z = (a.M1.x - a.M2.y) / 2.;
  if (a.dc >= 0) {
    a.l1r = (a.M1.x + a.M2.y) / 2. + sqrt(a.dc);
    a.l2r = (a.M1.x + a.M2.y) / 2. - sqrt(a.dc);
    a.l1i = 0;
    a.l2i = 0;
  } else {
    a.l1r = (a.M1.x + a.M2.y) / 2.;
    a.l2r = (a.M1.x + a.M2.y) / 2.;
    a.l1i = sqrt(-a.dc);
    a.l2i = -sqrt(-a.dc);
  }
  
  if (debug) {
    println("dc: " + a.dc);
    println("z: " + a.z);
  }
  
}
