class Affinity
{
  pt M1, M2, T; // Columns of a 2x3 matrix
  pt F; // fixed point
  
  Affinity() {
    M1 = new pt();
    M2 = new pt();
    T = new pt();
    F = new pt();
  }
  
  vec velocity(pt p) {
    return V(
      p.x * M1.x + p.y * M2.x + T.x,
      p.x * M1.y + p.y * M2.y + T.y);
  }
}
