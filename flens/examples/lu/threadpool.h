#ifndef LU_THREADPOOL_H
#define LU_THREADPOOL_H 1

#include <condition_variable>
#include <list>
#include <mutex>
#include <unordered_map>
#include <thread>

namespace flens {

class ThreadPool
{
    public:

        typedef std::function<void()>    Task;

        ThreadPool(int numThreads);

        void
        add(Task task);

        int
        join();

        void
        abort(int errorCode);

    private:

        Task
        getTask();

        void
        runThread();

        void
        taskDone();

        int                                         numThreads;
        int                                         numThreadsBusy;
        int                                         errorCode;

        std::list<Task>                           ready;

        std::condition_variable                     cv_master;
        std::condition_variable                     cv_thread;
        std::mutex                                  mtx;
};

} // namespace flens

#endif // LU_THREADPOOL_H
