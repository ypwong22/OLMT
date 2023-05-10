import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin

class MLPChained(BaseEstimator, RegressorMixin):
    def __init__(self, classifier, regressor):
        self.classifier = classifier
        self.regressor = regressor

    def fit(self, X, y):
      non_zero_indices = np.abs(y) > 1e-10
      self.classifier.fit(X, non_zero_indices)
      # non_zero_indices = np.where(self.classifier.predict(X) == 1)[0]
      if non_zero_indices.size > 0:
        self.regressor.fit(X[non_zero_indices], y[non_zero_indices])
      return self

    def predict(self, X):
      output = np.zeros(len(X))
      non_zero_indices = np.where(self.classifier.predict(X))[0]
      if non_zero_indices.size > 0:
          output[non_zero_indices] = self.regressor.predict(X[non_zero_indices])
      return output
