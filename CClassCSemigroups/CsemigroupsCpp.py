# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.1
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _CsemigroupsCpp
else:
    import _CsemigroupsCpp

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class SwigPyIterator(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _CsemigroupsCpp.delete_SwigPyIterator

    def value(self):
        return _CsemigroupsCpp.SwigPyIterator_value(self)

    def incr(self, n=1):
        return _CsemigroupsCpp.SwigPyIterator_incr(self, n)

    def decr(self, n=1):
        return _CsemigroupsCpp.SwigPyIterator_decr(self, n)

    def distance(self, x):
        return _CsemigroupsCpp.SwigPyIterator_distance(self, x)

    def equal(self, x):
        return _CsemigroupsCpp.SwigPyIterator_equal(self, x)

    def copy(self):
        return _CsemigroupsCpp.SwigPyIterator_copy(self)

    def next(self):
        return _CsemigroupsCpp.SwigPyIterator_next(self)

    def __next__(self):
        return _CsemigroupsCpp.SwigPyIterator___next__(self)

    def previous(self):
        return _CsemigroupsCpp.SwigPyIterator_previous(self)

    def advance(self, n):
        return _CsemigroupsCpp.SwigPyIterator_advance(self, n)

    def __eq__(self, x):
        return _CsemigroupsCpp.SwigPyIterator___eq__(self, x)

    def __ne__(self, x):
        return _CsemigroupsCpp.SwigPyIterator___ne__(self, x)

    def __iadd__(self, n):
        return _CsemigroupsCpp.SwigPyIterator___iadd__(self, n)

    def __isub__(self, n):
        return _CsemigroupsCpp.SwigPyIterator___isub__(self, n)

    def __add__(self, n):
        return _CsemigroupsCpp.SwigPyIterator___add__(self, n)

    def __sub__(self, *args):
        return _CsemigroupsCpp.SwigPyIterator___sub__(self, *args)
    def __iter__(self):
        return self

# Register SwigPyIterator in _CsemigroupsCpp:
_CsemigroupsCpp.SwigPyIterator_swigregister(SwigPyIterator)

class vector_long(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _CsemigroupsCpp.vector_long_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _CsemigroupsCpp.vector_long___nonzero__(self)

    def __bool__(self):
        return _CsemigroupsCpp.vector_long___bool__(self)

    def __len__(self):
        return _CsemigroupsCpp.vector_long___len__(self)

    def __getslice__(self, i, j):
        return _CsemigroupsCpp.vector_long___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _CsemigroupsCpp.vector_long___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _CsemigroupsCpp.vector_long___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _CsemigroupsCpp.vector_long___delitem__(self, *args)

    def __getitem__(self, *args):
        return _CsemigroupsCpp.vector_long___getitem__(self, *args)

    def __setitem__(self, *args):
        return _CsemigroupsCpp.vector_long___setitem__(self, *args)

    def pop(self):
        return _CsemigroupsCpp.vector_long_pop(self)

    def append(self, x):
        return _CsemigroupsCpp.vector_long_append(self, x)

    def empty(self):
        return _CsemigroupsCpp.vector_long_empty(self)

    def size(self):
        return _CsemigroupsCpp.vector_long_size(self)

    def swap(self, v):
        return _CsemigroupsCpp.vector_long_swap(self, v)

    def begin(self):
        return _CsemigroupsCpp.vector_long_begin(self)

    def end(self):
        return _CsemigroupsCpp.vector_long_end(self)

    def rbegin(self):
        return _CsemigroupsCpp.vector_long_rbegin(self)

    def rend(self):
        return _CsemigroupsCpp.vector_long_rend(self)

    def clear(self):
        return _CsemigroupsCpp.vector_long_clear(self)

    def get_allocator(self):
        return _CsemigroupsCpp.vector_long_get_allocator(self)

    def pop_back(self):
        return _CsemigroupsCpp.vector_long_pop_back(self)

    def erase(self, *args):
        return _CsemigroupsCpp.vector_long_erase(self, *args)

    def __init__(self, *args):
        _CsemigroupsCpp.vector_long_swiginit(self, _CsemigroupsCpp.new_vector_long(*args))

    def push_back(self, x):
        return _CsemigroupsCpp.vector_long_push_back(self, x)

    def front(self):
        return _CsemigroupsCpp.vector_long_front(self)

    def back(self):
        return _CsemigroupsCpp.vector_long_back(self)

    def assign(self, n, x):
        return _CsemigroupsCpp.vector_long_assign(self, n, x)

    def resize(self, *args):
        return _CsemigroupsCpp.vector_long_resize(self, *args)

    def insert(self, *args):
        return _CsemigroupsCpp.vector_long_insert(self, *args)

    def reserve(self, n):
        return _CsemigroupsCpp.vector_long_reserve(self, n)

    def capacity(self):
        return _CsemigroupsCpp.vector_long_capacity(self)
    __swig_destroy__ = _CsemigroupsCpp.delete_vector_long

# Register vector_long in _CsemigroupsCpp:
_CsemigroupsCpp.vector_long_swigregister(vector_long)

class vector_vector_long(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _CsemigroupsCpp.vector_vector_long_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _CsemigroupsCpp.vector_vector_long___nonzero__(self)

    def __bool__(self):
        return _CsemigroupsCpp.vector_vector_long___bool__(self)

    def __len__(self):
        return _CsemigroupsCpp.vector_vector_long___len__(self)

    def __getslice__(self, i, j):
        return _CsemigroupsCpp.vector_vector_long___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _CsemigroupsCpp.vector_vector_long___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _CsemigroupsCpp.vector_vector_long___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _CsemigroupsCpp.vector_vector_long___delitem__(self, *args)

    def __getitem__(self, *args):
        return _CsemigroupsCpp.vector_vector_long___getitem__(self, *args)

    def __setitem__(self, *args):
        return _CsemigroupsCpp.vector_vector_long___setitem__(self, *args)

    def pop(self):
        return _CsemigroupsCpp.vector_vector_long_pop(self)

    def append(self, x):
        return _CsemigroupsCpp.vector_vector_long_append(self, x)

    def empty(self):
        return _CsemigroupsCpp.vector_vector_long_empty(self)

    def size(self):
        return _CsemigroupsCpp.vector_vector_long_size(self)

    def swap(self, v):
        return _CsemigroupsCpp.vector_vector_long_swap(self, v)

    def begin(self):
        return _CsemigroupsCpp.vector_vector_long_begin(self)

    def end(self):
        return _CsemigroupsCpp.vector_vector_long_end(self)

    def rbegin(self):
        return _CsemigroupsCpp.vector_vector_long_rbegin(self)

    def rend(self):
        return _CsemigroupsCpp.vector_vector_long_rend(self)

    def clear(self):
        return _CsemigroupsCpp.vector_vector_long_clear(self)

    def get_allocator(self):
        return _CsemigroupsCpp.vector_vector_long_get_allocator(self)

    def pop_back(self):
        return _CsemigroupsCpp.vector_vector_long_pop_back(self)

    def erase(self, *args):
        return _CsemigroupsCpp.vector_vector_long_erase(self, *args)

    def __init__(self, *args):
        _CsemigroupsCpp.vector_vector_long_swiginit(self, _CsemigroupsCpp.new_vector_vector_long(*args))

    def push_back(self, x):
        return _CsemigroupsCpp.vector_vector_long_push_back(self, x)

    def front(self):
        return _CsemigroupsCpp.vector_vector_long_front(self)

    def back(self):
        return _CsemigroupsCpp.vector_vector_long_back(self)

    def assign(self, n, x):
        return _CsemigroupsCpp.vector_vector_long_assign(self, n, x)

    def resize(self, *args):
        return _CsemigroupsCpp.vector_vector_long_resize(self, *args)

    def insert(self, *args):
        return _CsemigroupsCpp.vector_vector_long_insert(self, *args)

    def reserve(self, n):
        return _CsemigroupsCpp.vector_vector_long_reserve(self, n)

    def capacity(self):
        return _CsemigroupsCpp.vector_vector_long_capacity(self)
    __swig_destroy__ = _CsemigroupsCpp.delete_vector_vector_long

# Register vector_vector_long in _CsemigroupsCpp:
_CsemigroupsCpp.vector_vector_long_swigregister(vector_vector_long)

class VecDouble(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _CsemigroupsCpp.VecDouble_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _CsemigroupsCpp.VecDouble___nonzero__(self)

    def __bool__(self):
        return _CsemigroupsCpp.VecDouble___bool__(self)

    def __len__(self):
        return _CsemigroupsCpp.VecDouble___len__(self)

    def __getslice__(self, i, j):
        return _CsemigroupsCpp.VecDouble___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _CsemigroupsCpp.VecDouble___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _CsemigroupsCpp.VecDouble___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _CsemigroupsCpp.VecDouble___delitem__(self, *args)

    def __getitem__(self, *args):
        return _CsemigroupsCpp.VecDouble___getitem__(self, *args)

    def __setitem__(self, *args):
        return _CsemigroupsCpp.VecDouble___setitem__(self, *args)

    def pop(self):
        return _CsemigroupsCpp.VecDouble_pop(self)

    def append(self, x):
        return _CsemigroupsCpp.VecDouble_append(self, x)

    def empty(self):
        return _CsemigroupsCpp.VecDouble_empty(self)

    def size(self):
        return _CsemigroupsCpp.VecDouble_size(self)

    def swap(self, v):
        return _CsemigroupsCpp.VecDouble_swap(self, v)

    def begin(self):
        return _CsemigroupsCpp.VecDouble_begin(self)

    def end(self):
        return _CsemigroupsCpp.VecDouble_end(self)

    def rbegin(self):
        return _CsemigroupsCpp.VecDouble_rbegin(self)

    def rend(self):
        return _CsemigroupsCpp.VecDouble_rend(self)

    def clear(self):
        return _CsemigroupsCpp.VecDouble_clear(self)

    def get_allocator(self):
        return _CsemigroupsCpp.VecDouble_get_allocator(self)

    def pop_back(self):
        return _CsemigroupsCpp.VecDouble_pop_back(self)

    def erase(self, *args):
        return _CsemigroupsCpp.VecDouble_erase(self, *args)

    def __init__(self, *args):
        _CsemigroupsCpp.VecDouble_swiginit(self, _CsemigroupsCpp.new_VecDouble(*args))

    def push_back(self, x):
        return _CsemigroupsCpp.VecDouble_push_back(self, x)

    def front(self):
        return _CsemigroupsCpp.VecDouble_front(self)

    def back(self):
        return _CsemigroupsCpp.VecDouble_back(self)

    def assign(self, n, x):
        return _CsemigroupsCpp.VecDouble_assign(self, n, x)

    def resize(self, *args):
        return _CsemigroupsCpp.VecDouble_resize(self, *args)

    def insert(self, *args):
        return _CsemigroupsCpp.VecDouble_insert(self, *args)

    def reserve(self, n):
        return _CsemigroupsCpp.VecDouble_reserve(self, n)

    def capacity(self):
        return _CsemigroupsCpp.VecDouble_capacity(self)
    __swig_destroy__ = _CsemigroupsCpp.delete_VecDouble

# Register VecDouble in _CsemigroupsCpp:
_CsemigroupsCpp.VecDouble_swigregister(VecDouble)

class VecVecdouble(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _CsemigroupsCpp.VecVecdouble_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _CsemigroupsCpp.VecVecdouble___nonzero__(self)

    def __bool__(self):
        return _CsemigroupsCpp.VecVecdouble___bool__(self)

    def __len__(self):
        return _CsemigroupsCpp.VecVecdouble___len__(self)

    def __getslice__(self, i, j):
        return _CsemigroupsCpp.VecVecdouble___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _CsemigroupsCpp.VecVecdouble___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _CsemigroupsCpp.VecVecdouble___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _CsemigroupsCpp.VecVecdouble___delitem__(self, *args)

    def __getitem__(self, *args):
        return _CsemigroupsCpp.VecVecdouble___getitem__(self, *args)

    def __setitem__(self, *args):
        return _CsemigroupsCpp.VecVecdouble___setitem__(self, *args)

    def pop(self):
        return _CsemigroupsCpp.VecVecdouble_pop(self)

    def append(self, x):
        return _CsemigroupsCpp.VecVecdouble_append(self, x)

    def empty(self):
        return _CsemigroupsCpp.VecVecdouble_empty(self)

    def size(self):
        return _CsemigroupsCpp.VecVecdouble_size(self)

    def swap(self, v):
        return _CsemigroupsCpp.VecVecdouble_swap(self, v)

    def begin(self):
        return _CsemigroupsCpp.VecVecdouble_begin(self)

    def end(self):
        return _CsemigroupsCpp.VecVecdouble_end(self)

    def rbegin(self):
        return _CsemigroupsCpp.VecVecdouble_rbegin(self)

    def rend(self):
        return _CsemigroupsCpp.VecVecdouble_rend(self)

    def clear(self):
        return _CsemigroupsCpp.VecVecdouble_clear(self)

    def get_allocator(self):
        return _CsemigroupsCpp.VecVecdouble_get_allocator(self)

    def pop_back(self):
        return _CsemigroupsCpp.VecVecdouble_pop_back(self)

    def erase(self, *args):
        return _CsemigroupsCpp.VecVecdouble_erase(self, *args)

    def __init__(self, *args):
        _CsemigroupsCpp.VecVecdouble_swiginit(self, _CsemigroupsCpp.new_VecVecdouble(*args))

    def push_back(self, x):
        return _CsemigroupsCpp.VecVecdouble_push_back(self, x)

    def front(self):
        return _CsemigroupsCpp.VecVecdouble_front(self)

    def back(self):
        return _CsemigroupsCpp.VecVecdouble_back(self)

    def assign(self, n, x):
        return _CsemigroupsCpp.VecVecdouble_assign(self, n, x)

    def resize(self, *args):
        return _CsemigroupsCpp.VecVecdouble_resize(self, *args)

    def insert(self, *args):
        return _CsemigroupsCpp.VecVecdouble_insert(self, *args)

    def reserve(self, n):
        return _CsemigroupsCpp.VecVecdouble_reserve(self, n)

    def capacity(self):
        return _CsemigroupsCpp.VecVecdouble_capacity(self)
    __swig_destroy__ = _CsemigroupsCpp.delete_VecVecdouble

# Register VecVecdouble in _CsemigroupsCpp:
_CsemigroupsCpp.VecVecdouble_swigregister(VecVecdouble)


def foo(a):
    return _CsemigroupsCpp.foo(a)

def foo2(gen):
    return _CsemigroupsCpp.foo2(gen)

def Pintar(v):
    return _CsemigroupsCpp.Pintar(v)

def gcd(a, b):
    return _CsemigroupsCpp.gcd(a, b)

def gcdL(v):
    return _CsemigroupsCpp.gcdL(v)

def prodEsc(*args):
    return _CsemigroupsCpp.prodEsc(*args)

def belongByGens(x, gen):
    return _CsemigroupsCpp.belongByGens(x, gen)

def computeMSG(generators):
    return _CsemigroupsCpp.computeMSG(generators)

def belongAxis(x, r):
    return _CsemigroupsCpp.belongAxis(x, r)

def axisIsSemigroup(gen, r):
    return _CsemigroupsCpp.axisIsSemigroup(gen, r)

def axisAreSemigroup(gen, setR):
    return _CsemigroupsCpp.axisAreSemigroup(gen, setR)

def MultiplicityAxis(generators, ray):
    return _CsemigroupsCpp.MultiplicityAxis(generators, ray)

def diamond(mult):
    return _CsemigroupsCpp.diamond(mult)

def pointBelongsDiamond(pt, eq):
    return _CsemigroupsCpp.pointBelongsDiamond(pt, eq)

def filterPoints(points, eq):
    return _CsemigroupsCpp.filterPoints(points, eq)

def eqRay(ray, hyperplanes):
    return _CsemigroupsCpp.eqRay(ray, hyperplanes)

def deleteRowZero(m):
    return _CsemigroupsCpp.deleteRowZero(m)

def affineTerm(eq, d):
    return _CsemigroupsCpp.affineTerm(eq, d)

def studyRays(rays, hyperplanes, integerDiamond, generators):
    return _CsemigroupsCpp.studyRays(rays, hyperplanes, integerDiamond, generators)

def existGenerator(equationsRay, affine, generators):
    return _CsemigroupsCpp.existGenerator(equationsRay, affine, generators)


