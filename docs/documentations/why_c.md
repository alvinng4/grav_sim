You may wonder why I chose C instead of Python. There are a few reasons:

1. Performance

    This one is obvious. As a low-level language, C is much faster than Python, and performance
    is essential for large-scale / long time scale simulations. Some may argue that Python could
    also be fast, with the help of Cython or Numba, but I just figured that it would be simpler
    and faster to just write in plain C. Besides, it is tempting to write fast but obscure Python
    code, which I have observed in some projects, but readability is a priority for me.

    At the very beginning, the project is written in Python with NumPy. As I rewrote the code
    in C, I have achieved a 400x - 1000x speedup. The difference has shrunk to 50x - 100x after 
    vectorizing the acceleration function to avoid the overhead from Python for loops.

2. Readability

    C is a small language with a small set of syntax. It doesn't have a lot of fancy
    features so it is easier to understand and maintain. In addition, everything is done
    explicitly in C, so there is no room for ambiguity.
    
3. Strong typing

    C is a strongly typed language and I found it easier to work with. In Python, I sometimes 
    still find bugs caused by wrong data types, even though I am already using type hints and
    mypy.

Nevertheless, Python is still useful in many ways. I would argue that it is better to use Python
when performance is not a concern. And it is really easy to plot graphs with Python.