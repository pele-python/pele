from __future__ import print_function
import numpy as np

class RunningMean(object):
    def __init__(self):
        self.count = 0
        self.sum = 0
    def __call__(self, value):
        try:
            self.sum += sum(value)
            self.count += len(value)
        except TypeError:
            self.sum += value
            self.count += 1
        return float(self.sum)/self.count
    
if __name__ == "__main__":
    mean = RunningMean()
    rand_list = np.random.randint(low=1, high=10, size=100)
    for i in rand_list:
        print(i, mean(i))
    add_list = np.random.randint(low=1, high=10, size=10)
    print(mean(add_list), float(sum(np.concatenate((rand_list, add_list))))/len(np.concatenate((rand_list, add_list))))
