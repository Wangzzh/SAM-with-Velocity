pts P = new pts();

void setup() {
  size(800, 800, P2D);
  P.declare();
  P.loadPts("data/pts");
}

void draw() {
  background(#FFFFFF);
  P.drawArrowsAndPoints();
}
