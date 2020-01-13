void compute_velocity_affinity(pt P1, pt P2, pt P3, vec V1, vec V2, vec V3, pt M1, pt M2, pt T) {
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
  M1.x = V1.x * i11 + V2.x * i21 + V3.x * i31;
  M2.x = V1.x * i12 + V2.x * i22 + V3.x * i32;
  T.x = V1.x * i13 + V2.x * i23 + V3.x * i33;
  M1.y = V1.y * i11 + V2.y * i21 + V3.y * i31;
  M2.y = V1.y * i12 + V2.y * i22 + V3.y * i32;
  T.y = V1.y * i13 + V2.y * i23 + V3.y * i33;
  
  if (debug) {
    println("P1: (" + P1.x + ", " + P1.y + ")");
    println("P2: (" + P2.x + ", " + P2.y + ")");
    println("P3: (" + P3.x + ", " + P3.y + ")");
    println("V1: (" + V1.x + ", " + V1.y + ")");
    println("V2: (" + V2.x + ", " + V2.y + ")");
    println("V3: (" + V3.x + ", " + V3.y + ")");
    println("Det: " + det);
    println("Inv:");
    println(i11 + " " + i12 + " " + i13);
    println(i21 + " " + i22 + " " + i23);
    println(i31 + " " + i32 + " " + i33);
    println("Affine:");
    println(M1.x + " " + M2.x + " " + T.x);
    println(M1.y + " " + M2.y + " " + T.y);
  }
}
