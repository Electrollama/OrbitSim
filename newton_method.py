from math import isnan

def zero(f, guess=0, tol=0, limit=None):
    """
    Return a guess of a zero of the function f, near guess.
    :param f: a differentiable function f(x)
    :param guess: initial guess for Newton's Method
    :param tol: tolerance, relatively close to zero
    :return: x such that |f(x)| <= tol
    """
    if not isinstance(guess, float):
        err_str = 'Incorrect type for optimization guess: {}'
        raise TypeError(err_str.format(type(guess)))
    return iter_improve(newton_update(f), lambda x: abs(f(x)) <= tol,
                        guess=guess, limit=limit)

def solution(f1, f2, guess=0, tol=0):
    """
    Find x such that f1(x) = f2(x)
    :param f1:
    :param f2:
    :param guess:
    :param tol:
    :return: x
    """
    return zero(lambda x: f1(x)-f2(x), guess, tol)

def iter_improve(update, done, guess=0.0, max_updates=1000, limit=None):
    """
    Iteratively improve guess with update until done returns a true value.
    :param update: function performing an update (ex. Newton's Method)
    :param done: boolean function returning True when the guess is sufficient
    :param guess: input for 'done'
    :param max_updates: maximum number of calls to 'update'
    :return: the first sufficient updated guess
    """
    count = 0
    while count < max_updates:
        old_guess = guess
        try:
            if done(guess):
                break # in the try statement in case a TypeError is thrown here
            guess = update(guess)
            if (guess is complex or
                isnan(guess)):
                raise TypeError
        except TypeError:
            warning_str = 'Value of type {} encountered during optimization.'
            warning_str = warning_str.format(type(old_guess))
            print(warning_str, guess)
            guess = old_guess*(1.0001)
        if limit is None:
            pass
        elif limit[1] is not None and guess > limit[1]:
            guess = limit[1]
        elif limit[0] is not None and guess < limit[0]:
            guess = limit[0]
        count += 1
    return guess

def approx_derivative(f, x, delta=1e-6):
    """
    Return an approximation to the derivative of f at x.
    :param f: function that is differentiable near x
    :param x: input to f
    :param delta: step size in x
    :return: df/dx
    """
    df = (f(x + delta) - f(x - delta))/2 #order 2 approximation
    return df/delta

def newton_update(f):
    """
    Return an update function for finding zeros using Newton's method.
    :param f: differentiable function
    :return: function to update guesses for x
    """
    return lambda x: x - f(x) / approx_derivative(f, x)
