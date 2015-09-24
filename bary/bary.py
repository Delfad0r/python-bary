import functools
import sympy


#exception
class BaryException(Exception):
	pass


#decorators
def _decorator_is_valid(func):
	@functools.wraps(func)
	def decorator(self, *args, **kwargs):
		if not self.is_valid():
			raise BaryException()
		return func(self, *args, **kwargs)
	return decorator


#enhanced tuple
class _EnhancedTuple:
	
	def __init__(self, *args, simplify = True):
		if simplify:
			self._tuple = tuple(sympy.simplify(ex) for ex in args)
			self.simplify()
		else:
			self._tuple = args
	
	def __iter__(self):
		return self._tuple.__iter__()
	
	def __len__(self):
		return self._tuple.__len__()
	
	def subs(self, *args, simplify = True, **kwargs):
		return type(self)(*(ex.subs(*args, **kwargs) for ex in self._tuple), simplify = simplify)
		
	def is_valid(self):
		return any(ex != 0 for ex in self._tuple)
	
	def is_vector(self):
		return sympy.Eq(sympy.simplify(sum(self._tuple)), 0) == True


#homogeneous tuple
class Htuple(_EnhancedTuple):
	
	def __repr__(self):
		return '(%s)' % ' : '.join(ex.__repr__() for ex in self._tuple)
	
	@_decorator_is_valid
	def simplify(self):
		self._tuple = tuple(sympy.factor(ex) for ex in self._tuple)
		m = functools.reduce(lambda x, y: x * sympy.denom(y), self._tuple, 1)
		self._tuple = tuple(ex * m for ex in self._tuple)
		d = functools.reduce(sympy.gcd, self._tuple)
		self._tuple = tuple(sympy.factor(ex / d) for ex in self._tuple)
	
	def is_normalized(self):
		return False

#normalized tuple
class Ntuple(_EnhancedTuple):
	
	def __init__(self, *args, simplify = True):
		super().__init__(*args, simplify = True)
	
	def __repr__(self):
		return '(%s)' % ' , '.join(ex.__repr__() for ex in self._tuple)
	
	@_decorator_is_valid
	def simplify(self):
		denom = sympy.simplify(sum(self._tuple))
		if denom == 0:
			raise BaryException()
		self._tuple = tuple(sympy.factor(sympy.simplify(ex / denom)) for ex in self._tuple)
	
	def is_normalized(self):
		return True
	

#class with many different transformations for points
class _PointTransform:
	
	def cycle(self):
		f, g, h = self._tuple
		subs_dict = {a : b, b : c, c : a}
		return type(self)(*(ex.subs(subs_dict, simultaneous = True, simplify = False) for ex in (h, f, g)))

	def swapxy(self):
		f, g, h = self._tuple
		subs_dict = {a : b, b : a}
		return type(self)(*(ex.subs(subs_dict, simultaneous = True, simplify = False) for ex in (g, f, h)))
		
	def swapyz(self):
		f, g, h = self._tuple
		subs_dict = {b : c, c : b}
		return type(self)(*(ex.subs(subs_dict, simultaneous = True, simplify = False) for ex in (f, h, g)))

	def swapzx(self):
		f, g, h = self._tuple
		subs_dict = {c : a, a : c}
		return type(self)(*(ex.subs(subs_dict, simultaneous = True, simplify = False) for ex in (h, g, f)))
	
	def homotetic(self, p, h):
		if not p.is_normalized():
			p = Npoint(*p)
		me = self
		if not me.is_normalized():
			me = Npoint(*me)
		lp, lme = tuple(p), tuple(me)
		return type(self)(*(h * (lme[i] - lp[i]) + lp[i] for i in (0, 1, 2)))
	
	def isotomic(self):
		return type(self)(*(1 / ex for ex in self._tuple))

	def isogonal(self):
		lp = self._tuple
		return type(self)(a ** 2 / lp[0], b ** 2 / lp[1], c ** 2 / lp[2])


#homogeneous point
class Hpoint(Htuple, _PointTransform):
	pass

#normalized point
class Npoint(Ntuple, _PointTransform):
	pass

#homogeneous vector
class Hvector(Htuple, _PointTransform):
	pass

#normalized vector
class Nvector(Ntuple, _PointTransform):
	
	def simplify(self):
		self._tuple = tuple(sympy.factor(sympy.simplify(ex)) for ex in self._tuple)


#curve described by an homogeneous equation
class Curve:
	
	def __init__(self, eqn, simplify = True):
		if isinstance(eqn, sympy.Equality):
			self._eqn = eqn.lhs - eqn.rhs
		else:
			self._eqn = eqn
		if simplify:
			self._eqn = sympy.simplify(self._eqn)
			self.simplify()
	
	def __repr__(self):
		self.simplify()
		return self.eqn().__repr__()
	
	def eqn(self):
		return sympy.Eq(self._eqn, 0)
	
	def subs(self, *args, simplify = True, **kwargs):
		return Curve(self._eqn.subs(*args, **kwargs), simplify = simplify)
	
	def cycle(self):
		subs_dict = {a : b, b : c, c : a, x : y, y : z, z : x}
		return self.subs(subs_dict, simultaneous = True, simplify = False)

	def swapxy(self):
		subs_dict = {a : b, b : a, x : y, y : x}
		return self.subs(subs_dict, simultaneous = True, simplify = False)
		
	def swapyz(self):
		subs_dict = {b : c, c : b, y : z, z : y}
		return self.subs(subs_dict, simultaneous = True, simplify = False)

	def swapzx(self):
		subs_dict = {c : a, a : c, z : x, x : z}
		return self.subs(subs_dict, simultaneous = True, simplify = False)
	
	def homotetic(self, p, h):
		s = x + y + z
		if not p.is_normalized():
			p = Npoint(*p)
		lp = tuple(p)
		subs_dict = {x : h * (x / s - lp[0]) + lp[0], y : h * (y / s - lp[1]) + lp[1], z : h * (z / s - lp[2]) + lp[2]}
		return self.subs(subs_dict, simultaneous = True)
	
	def isotomic(self):
		subs_dict = {x : 1 / x, y : 1 / y, z : 1 / z}
		return self.subs(subs_dict, simultaneous = True)
	
	def isogonal(self):
		subs_dict = {x : a ** 2 / x, y : b ** 2 / y, z : c ** 2 / z}
		return self.subs(subs_dict, simultaneous = True)
	
	def is_valid(self):
		return self._eqn not in (True, False)
	
	@_decorator_is_valid
	def simplify(self):
		self._eqn = sympy.factor(sympy.numer(sympy.factor(self._eqn)))
		if isinstance(self._eqn, sympy.Mul):
			def is_not_constant(ex):
				sym = ex.free_symbols
				return x in sym or y in sym or z in sym
			self._eqn = sympy.Mul(*(filter(is_not_constant, self._eqn.args)))


#homogeneous solve
def hsolve(f, *symbols, **flags):
	flags['dict'] = True
	solutions = sympy.solve(f, *symbols, **flags)
	ret = []
	for sol in solutions:
		freedom_deg = 0
		for sym in symbols:
			if sym not in sol:
				freedom_deg += 1
				sol[sym] = sym
			else:
				sol[sym] = sol[sym].subs([(sympy.Abs(s), s) for s in symbols])
		if freedom_deg != 1:
			raise BaryException()
		ret.append(Htuple(*tuple(sol[sym] for sym in symbols)))
	return ret


#utility
def to_list_of_equations(*args):
	ret = []
	for ex in args:
		if isinstance(ex, (list, tuple)):
			ret += to_list_of_equations(*(x for x in ex))
		elif isinstance(ex, sympy.And):
			ret += to_list_of_equations(*ex.args)
		else:
			ret.append(ex)
	return ret

def cyclic3uple(ex0):
	subs_dict = {a : b, b : c, c : a}
	ex1 = ex0.subs(subs_dict, simultaneous = True)
	ex2 = ex1.subs(subs_dict, simultaneous = True)
	return (ex0, ex1, ex2)

def cycle(obj):
	return obj.cycle()

def swapxy(o):
	return obj.swapxy()
	
def swapyz(o):
	return obj.swapyz()

def swapzx(o):
	return obj.swapzx()

def homothetic(p, h, obj):
	return obj.homotetic(p, h)

def isotomic(obj):
	return obj.isotomic()

def isogonal(obj):
	return obj.isogonal()


#points
def Point(x, y, z):
	return Hpoint(x, y, z)

def is_same_point(p1, p2):
	l1, l2 = tuple(p1), tuple(p2)
	return sympy.simplify(sympy.And(sympy.Eq(l1[0] * l2[1] - l1[1] * l2[0], 0), sympy.Eq(l1[0] * l2[2] - l1[2] * l2[0])))

def collinear(p1, p2, p3):
	l1, l2, l3 = list(p1), list(p2), list(p3)
	return sympy.simplify(sympy.Eq(sympy.det(sympy.Matrix([l1, l2, l3])), 0))

def distance(p1, p2):
	return vector_length(displacement_vector(p1, p2, normalized = True))

def area(p1, p2, p3):
	if not p1.is_normalized():
		p1 = Npoint(*p1)
	if not p2.is_normalized():
		p2 = Npoint(*p2)
	if not p3.is_normalized():
		p3 = Npoint(*p3)
	return S / 2 * sympy.det(sympy.Matrix([list(p1), list(p2), list(p3)]))

def midpoint(p1, p2, weight1 = 1, weight2 = 1):
	if not p1.is_normalized():
		p1 = Npoint(*p1)
	if not p2.is_normalized():
		p2 = Npoint(*p2)
	l1, l2 = tuple(p1), tuple(p2)
	return Hpoint(*(l1[i] * weight1 + l2[i] * weight2 for i in (0, 1, 2)))


#lines
def Line(l, m, n):
	return Curve(l * x + m * y + n * z)

def line_through_two_points(p1, p2):
	l1, l2 = list(p1), list(p2)
	return Curve(sympy.det(sympy.Matrix([l1, l2, [x, y, z]])))

def infinity_point(r):
	p = intersect(r, Line(1, 1, 1))
	if len(p) != 1:
		raise BaryException()
	return p[0]

def perpendicular_through_point(r, p):
	u, v, w = tuple(sympy.Dummy(s) for s in ('u', 'v', 'w'))
	sol = hsolve([u + v + w, vector_perpendicular(Hpoint(u, v, w), infinity_point(r))], u, v, w)
	if len(sol) != 1:
		raise BaryException()
	return line_through_two_points(p, sol[0])

def perpendicular(r1, r2):
	return vector_perpendicular(infinity_point(r1),	infinity_point(r2))

def parallel_through_point(r, p):
	return line_through_two_points(p, infinity_point(r))

def parallel(r1, r2):
	return concur(r1, r2, Line(1, 1, 1))
	

#circles
def Circle(u, v, w):
	return Curve(-a ** 2 * y * z - b ** 2 * z * x - c ** 2 * x * y + (x + y + z) * (u * x + v * y + w * z))

def circle_through_three_points(p1, p2, p3):
	u, v, w = tuple(sympy.Dummy(s) for s in ('u', 'v', 'w'))
	circle = Circle(u, v, w)
	equations = [belongs_to(circle, p) for p in (p1, p2, p3)]
	return circle.subs(sympy.solve(equations, u, v, w))

def circle_center_radius(c, r):
	return Curve(distance(c, Point(x, y, z)) ** 2 - r ** 2)

def pow(circle, p):
	if not p.is_normalized():
		p = Npoint(*p)
	l = tuple(p)
	return circle.eqn().lhs.subs([(x, l[0]), (y, l[1]), (z, l[2])])
	
	
#curves
def intersect(c1, c2):
	return [Point(*sol) for sol in hsolve([c1.eqn(), c2.eqn()], x, y, z)]

def concur(c1, c2, c3):
	points = intersect(c1, c2)
	if len(points) == 0:
		return False
	return sympy.Or(*(belongs_to(c3, p) for p in points))

def belongs_to(c, p):
	l = tuple(p)
	return sympy.simplify(c.eqn().subs([(x, l[0]), (y, l[1]), (z, l[2])]))


#vectors
def displacement_vector(p1, p2, normalized = False):
	if not p1.is_normalized():
		p1 = Npoint(*p1)
	if not p2.is_normalized():
		p2 = Npoint(*p2)
	l1, l2 = tuple(p1), tuple(p2)
	v = (l2[i] - l1[i] for i in (0, 1, 2))
	if normalized:
		return Nvector(*v)
	else:
		return Hvector(*v)

def vector_length(v):
	x0, y0, z0 = tuple(v)
	return sympy.simplify(sympy.sqrt(-(a ** 2 * y0 * z0 + b ** 2 * z0 * x0 + c ** 2 * x0 * y0)))

def vector_perpendicular(v1, v2):
	l1, l2 = tuple(v1), tuple(v2)
	return sympy.Eq(sympy.simplify(
		a ** 2 * (l1[1] * l2[2] + l1[2] * l2[1]) +
		b ** 2 * (l1[2] * l2[0] + l1[0] * l2[2]) +
		c ** 2 * (l1[0] * l2[1] + l1[1] * l2[0]), 0))


#symbols

#coordinates
x, y, z = sympy.symbols('x y z', real = True)

#sides of the triangle (> 0)
a, b, c = sympy.symbols('a b c', positive = True)

#surface-related
S_A = (-a * a + b * b + c * c) / 2
S_B = ( a * a - b * b + c * c) / 2
S_C = ( a * a + b * b - c * c) / 2
S_AB = S_A * S_B
S_BC = S_B * S_C
S_CA = S_C * S_A
S = sympy.sqrt(S_AB + S_BC + S_CA)

