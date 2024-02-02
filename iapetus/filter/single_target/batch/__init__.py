"""Batch least squares module."""


class BatchIterator:
    """Batch processing iterator class."""

    def __init__(self, iter_tol: float, max_iter: int):
        """
        Args:
            iter_tol (float): Iteration error tolerance; determines when
                batch processor iterations have converged
            max_iter (int): Maximum number times to iterate on batch processor
                solution; guards agains endless loops in presence of
                non-convergence
        """
        self.iter_tol = iter_tol
        self.max_iter = max_iter
