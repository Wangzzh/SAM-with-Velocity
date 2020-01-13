boolean debug = false;

boolean showField = false;

pts P = new pts();
Affinity A = new Affinity();

void setup() {
  size(800, 800, P2D);
  P.declare();
  P.loadPts("data/pts");
}

void draw() {
  background(#FFFFFF);
  if (P.nv >= 6) {
    compute_velocity_affinity(P.G[0], P.G[2], P.G[4], V(P.G[0], P.G[1]), V(P.G[2], P.G[3]), V(P.G[4], P.G[5]), A.M1, A.M2, A.T);
  }
  stroke(black); strokeWeight(3);
  P.drawArrowsAndPoints();  
  
  if (showField) {
    stroke(blue); strokeWeight(1);
    for (float i = 0; i <= width; i += 50) {
      for (float j = 0; j <= height; j += 50) {
        pt p = P(i, j);
        vec v = A.velocity(p);
        arrow(p, v);
      }
    }
  }
}
