import numpy as np

def sign_pph(y):
    """

    Parameters
    ----------
    y : Value.

    Returns
    -------
    1.0 if np.sign(y) > 0
    -1.0 else.

    """
    y = np.array(y)
    if np.sign(y) > 0:
        return 1.0;
    else:
        return -1.0;