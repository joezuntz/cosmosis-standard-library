import numpy as np
import timeit

def read_uncompressed(filename):
    """
    Read the standard uncompressed form of this
    covariance matrix, which has the size of the
    matrix as an integer on the first line and then
    a line for each element, including the upper and
    lower diagonals.

    Parameters
    ----------
    filename: str

    Returns
    -------
    C: 2d array
        The dense form of the matrix
    """
    cov_data = np.loadtxt(filename)
    n = int(cov_data[0])
    C = cov_data[1:].reshape((n, n))
    return C


def read_compressed(filename):
    """
    Read the compressed form of the matrix made by
    the compress function below, which stores only
    the upper triangle.

    Parameters
    ----------
    filename: str

    Returns
    -------
    C: 2d array
        The dense form of the matrix
    """
    cov_data = np.loadtxt(filename)
    n = int(cov_data[0])
    C = np.zeros((n, n))

    # Fill in the upper triangle
    ij = np.triu_indices(n)
    C[ij] = cov_data[1:]

    # Flip and add the matrix to fill in the
    # lower triangle, and then subtract off
    # the diagonal to stop it being added twice.
    C = C + C.T - np.diag(C.diagonal())
    return C


def compress(filename_in, filename_out):
    """
    Compress the input format
    """
    C = read_uncompressed(filename_in)
    n = C.shape[0]

    # indices for the upper triangle of the matrix
    ij = np.triu_indices(n)

    # Save the size as well just to be slightly consistent with
    # the input format
    data = np.concatenate([[n], C[ij]])
    np.savetxt(filename_out, data)

def main():
    filename_in = "Pantheon+SH0ES_STAT+SYS.cov"
    filename_out = "Pantheon+SH0ES_STAT+SYS.cov_compressed.gz"
    compress(filename_in, filename_out)

    # Check that the compression worked
    t1 = timeit.default_timer()
    C1 = read_uncompressed(filename_in)
    t2 = timeit.default_timer()
    C2 = read_compressed(filename_out)
    t3 = timeit.default_timer()
    print(C1 - C2)
    print("Matrices equal.")
    print("Time to read uncompressed:", t3 - t2)
    print("Time to read compressed:", t2 - t1)

if __name__ == '__main__':
    main()