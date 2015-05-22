#include <flens/examples/lu/threadpool.h>
#include <cassert>
#include <cstdlib>
#include <fstream>

namespace flens {

ThreadPool::ThreadPool(int numThreads_)
    : numThreads(numThreads_), numThreadsBusy(0), errorCode(0)
{
    for (int i=0; i<numThreads; ++i) {
        auto t = std::thread( [=] { this->runThread(); } );
        t.detach();
    }
}

void
ThreadPool::add(Task task)
{
    {
        std::unique_lock<std::mutex>  lock(mtx);

        if (errorCode!=0) {
            return;
        }
        ready.push_back(task);
    }
    cv_thread.notify_one();
}

void
ThreadPool::abort(int errorCode_)
{
    assert(errorCode_!=0);
    {
        std::unique_lock<std::mutex> lock(mtx);

        errorCode = errorCode_;
    }
}

int
ThreadPool::join()
{
    std::unique_lock<std::mutex> lock(mtx);

    while (ready.size()+numThreadsBusy>0) {
        cv_master.wait(lock);
    }

    // clean up and restore error code
    ready.clear();
    int errorCode_ = errorCode;
    errorCode = 0;

    return errorCode_;
}

ThreadPool::Task
ThreadPool::getTask()
{
    Task task;
    bool jobReady = false;
    {
        std::unique_lock<std::mutex> lock(mtx);

        while (ready.size()==0) {
            cv_thread.wait(lock);
        }

        task = ready.front();
        ready.pop_front();
        ++numThreadsBusy;

        if (ready.size()>0) {
            jobReady = true;
        }
    }
    if (jobReady) {
        cv_thread.notify_one();
    }
    return task;
}

void
ThreadPool::runThread()
{
    while(1) {
        Task task = getTask();
        task();
        taskDone();
    }
}

void
ThreadPool::taskDone()
{
    bool allDone = false;

    {
        std::unique_lock<std::mutex> lock(mtx);

        if (errorCode!=0) {
            allDone= true;
        }
        --numThreadsBusy;
        if (ready.size()==0) {
            allDone = true;
        }
    }
    if (allDone) {
        cv_master.notify_one();
    }
}

} // namespace flens
