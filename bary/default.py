from bary.bary import *

#important points
A = Point(1, 0, 0)
B = Point(0, 1, 0)
C = Point(0, 0, 1)

G = Point(1, 1, 1)
I = Point(a, b, c)
H = Point(S_BC, S_CA, S_AB)
O = isogonal(H)
I_a = Point(-a, b, c)
I_b = Point(a, -b, c)
I_c = Point(a, b, -c)

#important lines
AB = BA = Line(0, 0, 1)
BC = CB = Line(1, 0, 0)
CA = AC = Line(0, 1, 0)
infinity_line = Line(1, 1, 1)
euler_line = line_through_two_points(O, G)

#important circles
circumcircle = Circle(0, 0, 0)
incircle = Circle(*cyclic3uple(((-a + b + c) / 2) ** 2))
excircle_A = Circle(((a + b + c) / 2) ** 2, ((a + b - c) / 2) ** 2, ((a - b + c) / 2) ** 2)
excircle_B = excircle_A.cycle()
excircle_C = excircle_B.cycle()
