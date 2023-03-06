# python -m doctest param.py

class param(float):
    '''
    Parameter with uncertainty property.

    Examples
    --------
    >>> p = param(5, 0)
    >>> p
    5 ± 0
    >>> f'{p:06.3f}'
    '05.000'
    >>> p + 2
    7.0
    '''
    def __new__(cls, value, unc=None):
        instance = super().__new__(cls, value)
        instance.unc = unc
        instance.value = value
        return instance

    def __repr__(self):
        return f'{self.value}' + ('' if self.unc is None  else f' ± {self.unc}')


class nesteddict(dict):
    '''
    A named and nested dictionary with multi-dimensional indexing.

    Example
    -------
    >>> d = nesteddict(norm=[1, 2, 3])
    >>> d.rv = 1.5
    >>> d.wave = [1, 3., 4.]
    >>> d.bkg = param(0.9, 0)
    >>> d.atm = {'H2O': 0.9, 'O2': 0.8, 'rv': 5}
    >>> d
    norm: [1, 2, 3]
    rv: 1.5
    wave: [1, 3.0, 4.0]
    bkg: 0.9 ± 0
    atm: {'H2O': 0.9, 'O2': 0.8, 'rv': 5}

    >>> d.flat()
    {('norm', 0): 1, ('norm', 1): 2, ('norm', 2): 3, 'rv': 1.5, ('wave', 0): 1, ('wave', 1): 3.0, ('wave', 2): 4.0, 'bkg': 0.9 ± 0, ('atm', 'H2O'): 0.9, ('atm', 'O2'): 0.8, ('atm', 'rv'): 5}

    >>> d['wave', 0]
    1

    >>> d['wave', 0] = 5
    '''
    __getattr__ = dict.__getitem__
    values = lambda _: [*super().values()]   # return a simple list instead of dict_values object
    keys = lambda _: [*super().keys()]   # return a simple list instead of dict_values object

    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __getitem__(self, key):
        # allow tuple indexing
        if isinstance(key, tuple):
            return self[key[0]][key[1]]
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            self[key[0]][key[1]] = value
        else:
            super().__setitem__(key, value)

    def __setattr__(self, key, value):
        self[key] = value

    def update(self, *args, **kwargs):
        # force calling _as_param via setitem
        for d in args + (kwargs,):
            for k in d:
                self[k] = d[k]

    def flat(self):
        d = {}
        for key, values in self.items():
            if isinstance(values, list):
                # numerate names
                for (i, val) in enumerate(values):
                    d[(key, i)] = val
            elif isinstance(values, dict):
                for (k, val) in values.items():
                    d[(key, k)] = val
            else:
               d[key] = values
        return d

    def __add__(self, d):
        p = self.__class__(self, d)
        return p

    def __repr__(self):
        return "\n".join([f'{k}: '+repr(v).replace('\n',', ') for k,v in self.items()])


class Params(nesteddict):
    '''
    A collection/group of param

    Example
    -------
    >>> par = Params(Params(rv=(1.5, 0)), norm=[1, 2, 3])
    >>> par.wave = [1, 3., 4.]
    >>> par.bkg = param(0.9, 0)
    >>> par.atm = {'H2O': 0.9, 'O2': param(0.8, 0.3), 'rv': 5}
    >>> par['wave', 0]
    1

    >>> par['wave', 0] = 5
    >>> par['atm', 'O2'].unc = 0.3

    >>> par.flat()
    {'rv': 1.5 ± 0, ('norm', 0): 1, ('norm', 1): 2, ('norm', 2): 3, ('wave', 0): 5, ('wave', 1): 3.0, ('wave', 2): 4.0, 'bkg': 0.9 ± 0, ('atm', 'H2O'): 0.9, ('atm', 'O2'): 0.8 ± 0.3, ('atm', 'rv'): 5}

    >>> par.update({'b': 5})
    >>> par.vary().keys()
    dict_keys([('norm', 0), ('norm', 1), ('norm', 2), ('wave', 0), ('wave', 1), ('wave', 2), ('atm', 'H2O'), ('atm', 'O2'), ('atm', 'rv'), 'b'])
    >>> vary = ['rv', ('norm', 0), ('wave', 1)]
    >>> vals = [901, 902, 903]
    >>> par.update(dict(zip(vary, vals)))
    >>> par
    rv: 901
    norm: [902, 2, 3]
    wave: [5, 903, 4.0]
    bkg: 0.9 ± 0
    atm: H2O: 0.9, O2: 0.8 ± 0.3, rv: 5
    b: 5

    >>> parcopy = Params(par)

    >>> par.b = [(6., 0)]
    >>> par + {('norm', 2): -1, 'c': 77}
    rv: 901
    norm: [902, 2, -1]
    wave: [5, 903, 4.0]
    bkg: 0.9 ± 0
    atm: H2O: 0.9, O2: 0.8 ± 0.3, rv: 5
    b: [6.0 ± 0]
    c: 77
    '''
    def __setitem__(self, key, value):
        super().__setitem__(key, self._as_param(value))

    def _as_param(self, value):
        if isinstance(value, (param, Params)):
             p = value
        elif isinstance(value, (float, int)):
             p = param(value)
        elif isinstance(value, tuple):
             p = param(*value)
        elif type(value).__name__ in ('list', 'ndarray'):
             p = [self._as_param(val) for val in value]
        elif isinstance(value, dict):
             p = Params(value)
        else:
            print(value, type(value), 'not supported')
        return p

    def vary(self):
        return {k: v for k,v in self.flat().items() if v.unc is not 0}
