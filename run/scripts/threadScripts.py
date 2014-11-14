#!/usr/bin/python
#Author: Lori Glenwinkel
import Queue
import threading
import os

queue = Queue.Queue()

class ThreadScripts(threading.Thread):
    """Threaded script execution"""
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue

    def run(self):
        while True:
            #grabs command from queue
            command=self.queue.get()
            #execute program
            os.system(command)
            #signal to queue job is done
            self.queue.task_done()

def spawn_threads(NumOfThreads):
    """#spawn a pool of threads, and pass them queue instance """
    for i in range(NumOfThreads):
        t = ThreadScripts(queue)
        t.setDaemon(True)
        t.start()
    
def pop_queue(data_queue):
    """populate queue with data"""
    for n in data_queue:
        queue.put(n)

def execute_queue(queueN):
    """send items in queue to be executed on different threads in parallel"""
    spawn_threads(len(queueN))
    pop_queue(queueN)
    queue.join()
