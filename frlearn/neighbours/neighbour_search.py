"""Nearest neighbour searches"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Callable, Union

import numpy as np
from sklearn.neighbors._unsupervised import NearestNeighbors


class NNSearch(ABC):
    """
    Abstract base class for nearest neighbour searches. Subclasses must
    implement __init__, construct and Index.query.
    """

    @abstractmethod
    def __init__(self, **kwargs):
        pass

    @abstractmethod
    def construct(self, X) -> Index:
        """
        Construct the index based on the data X.

        Parameters
        ----------
        X : array shape=(n_instances, n_features, )
            Construction instances.

        Returns
        -------
        I : Index
            Constructed index
        """
        index = self.Index.__new__(self.Index)
        index._X = X
        index._len = len(X)
        return index

    class Index(ABC):

        _X: np.array
        _len: int

        def query_self(self, k: Union[int, float, None]):
            if callable(k):
                k = k(len(self) - 1)
            return [a[:, 1:] for a in self.query(self._X, k + 1)]

        def query(self, X, k: Union[int, Callable[[int], int]]):
            """
            Identify the k nearest neighbours for each of the instances in X.

            Parameters
            ----------
            X : array shape=(n_instances, n_features, )
                Query instances.

            k : int or (int -> int)
                Number of neighbours to return. Should be either a positive integer not larger than the index size,
                or a function that takes the size of the index and returns such an integer.

            Returns
            -------
            I : array shape=(n_instances, k, )
                Indices of the k nearest neighbours among the construction
                instances for each query instance.

            D : array shape=(n_instances, k, )
                Distances to the k nearest neighbours among the construction
                instances for each query instance.
            """
            if callable(k):
                k = k(len(self))
            return self._query(X, k)

        @abstractmethod
        def _query(self, X, k: int):
            pass

        def __len__(self):
            return self._len


class BallTree(NNSearch):
    """
    Nearest neighbour search with a Ball tree.

    Parameters
    ----------
    metric : str, default='manhattan'
        The metric through which distances are defined.

    leaf_size : int, default=30
        The leaf size to be used for the Ball tree.

    n_jobs : int, default=1
        The number of parallel jobs to run for neighbour search. -1 means using
        all processors.
    """

    def __init__(self, *, metric: str = 'manhattan', leaf_size: int = 30,
                 n_jobs: int = 1):
        self.construction_params = {
            'algorithm': 'ball_tree',
            'metric': metric,
            'leaf_size': leaf_size,
            'n_jobs': n_jobs,
        }

    def construct(self, X) -> Index:
        index = super().construct(X)
        index.tree = NearestNeighbors(**self.construction_params).fit(X)
        return index

    class Index(NNSearch.Index):

        tree: NearestNeighbors

        def _query(self, X, k: int):
            return self.tree.kneighbors(X, n_neighbors=k)[::-1]


class KDTree(NNSearch):
    """
    Nearest neighbour search with a KD-tree.

    Parameters
    ----------
    metric : str, default='manhattan'
        The metric through which distances are defined.

    leaf_size : int, default=30
        The leaf size to be used for the KD-tree.

    n_jobs : int, default=1
        The number of parallel jobs to run for neighbour search. -1 means using
        all processors.
    """

    def __init__(self, *, metric: str = 'manhattan', leaf_size: int = 30,
                 n_jobs: int = 1):
        self.construction_params = {
            'algorithm': 'kd_tree',
            'metric': metric,
            'leaf_size': leaf_size,
            'n_jobs': n_jobs,
        }

    def construct(self, X) -> Index:
        index = super().construct(X)
        index.tree = NearestNeighbors(**self.construction_params).fit(X)
        return index

    class Index(NNSearch.Index):

        tree: NearestNeighbors

        def _query(self, X, k: int):
            return self.tree.kneighbors(X, n_neighbors=k)[::-1]
