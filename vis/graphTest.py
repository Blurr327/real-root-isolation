from manimlib import *

def readPoly(pf):
    poly = []
    with open(pf) as f:
    
        for l in f:
            poly.append(int(l.rstrip()))
    return poly
p1 = readPoly("poly")
p2 = readPoly("poly2")


def evalPoly(x, p):
    res = 0
    for i in range(len(p)):
        res += (pow(x,i) * p[i])
    return res

class PolyZoomScene(Scene):
    def construct(self):
        #p1 = [0, 1, 0, 0.1]   # small scale
        #p2 = [5, -2, 0.5]     # larger scale

        # One fixed axes
        axes = Axes(
            x_range=(-50, 50, 1),
            y_range=(-100, 100, 2),
            axis_config={
                "include_numbers": True,
            }
        )

        for number in axes.x_axis.numbers:
            number.scale(2.0)

        for number in axes.y_axis.numbers:
            number.scale(2.0)

        graph1 = axes.get_graph(lambda x: evalPoly(x, p1), color=BLUE)
        graph2 = axes.get_graph(lambda x: evalPoly(x, p2), color=BLUE)

        self.play(ShowCreation(axes))
        self.play(ShowCreation(graph1))
        self.wait()

        # --- Camera zoom OUT to fit bigger polynomial ---
        frame = self.camera.frame
        
        scale_ratio = 3 #change this to correct how much the graph gets scales between each poly
        self.play(
            frame.animate.scale(scale_ratio).move_to(axes.get_center()),
            run_time=2
        )

        # --- Transition graph ---
        self.play(
            ReplacementTransform(graph1, graph2),
            run_time=2
        )

        self.wait()

