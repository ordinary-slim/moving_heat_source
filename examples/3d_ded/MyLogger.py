import time

class IterationData:
    def __init__(self):
        self.isCoupled  = False
        self.startTime = 0.0
        self.endTime = +1e9
        self.duration = +1e9

class MyLogger:
    def __init__(self):
        self.iterations = []
    def iterate(self, driver):
        iterData = IterationData()
        iterData.startTime = time.time()
        driver.iterate()
        iterData.endTime = time.time()
        iterData.duration = iterData.endTime - iterData.startTime
        iterData.isCoupled = driver.isCoupled
        self.iterations.append( iterData )
