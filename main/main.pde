boolean debug = false;

boolean showField = false;
boolean showTrail = true;
boolean showVelocity = true;
boolean animate = false;
float t = 1, dt = 0.01;

pts P = new pts();
pts Q = new pts();
Affinity A = new Affinity();

void setup() {
  size(800, 800, P2D);
  P.declare();
  P.loadPts("data/pts");
  Q.declare();
  Q.nv = 3;
}

void draw() {
  background(#FFFFFF);
  if (animate) {
    t += dt;
    if (t >= 1.0) t = 0;
  }
  
  if (P.nv >= 7) {
    compute_velocity_affinity(P.G[0], P.G[2], P.G[4], V(P.G[0], P.G[1]), V(P.G[2], P.G[3]), V(P.G[4], P.G[5]), A);
    stroke(red); strokeWeight(3); fill(red);
    ellipse(A.F.x, A.F.y, 10, 10);
    
    stroke(red); strokeWeight(3); noFill();
    Q.G[0] = A.apply(P.G[0], t);
    Q.G[1] = A.apply(P.G[2], t);
    Q.G[2] = A.apply(P.G[4], t);
    if (debug) {
      println("Q1: (" + Q.G[0].x + ", " + Q.G[0].y + ")");
      println("Q2: (" + Q.G[1].x + ", " + Q.G[1].y + ")");
      println("Q3: (" + Q.G[2].x + ", " + Q.G[2].y + ")");
    }
    show(Q.G[0], Q.G[1], Q.G[2]);
    
    stroke(black); strokeWeight(3); noFill();
    show(P.G[0], P.G[2], P.G[4]);
    
  }
  stroke(black); strokeWeight(3);
  P.drawArrowsAndPoints();  
  
  if (showField && P.nv >= 7) {
    stroke(blue); strokeWeight(1);
    for (float i = 0; i <= width; i += 50) {
      for (float j = 0; j <= height; j += 50) {
        pt p = P(i, j);
        vec v = A.velocity(p);
        arrow(p, v);
      }
    }
  }
  
  if (showVelocity && P.nv >= 7) {
    stroke(dgreen); strokeWeight(3);
    vec v = A.velocity(P.G[6]);
    arrow(P.G[6], v);
  }
  
  if (showTrail && P.nv >= 7) {
    stroke(orange); strokeWeight(3);
    beginShape(); 
    for (float x = 0; x <= 3.0; x += 0.02) {
      pt p1 = A.apply(P.G[6], x);
      vertex(p1.x, p1.y);
    }
    endShape(OPEN);
  }
}
