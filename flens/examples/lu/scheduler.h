#ifndef LU_SCHEDULER_H
#define LU_SCHEDULER_H 1

#include <condition_variable>
#include <list>
#include <unordered_map>
#include <mutex>
#include <thread>
#include <vector>

#include <flens/examples/lu/format.h>

namespace flens {

class Scheduler
{
    public:

        ///
        ///  Internally task are identified by an integer.
        ///
        typedef long                     Key;
        ///
        ///  And a task is just a function pointer.
        ///
        typedef std::function<void()>    Task;

        ///
        ///  The scheduler and thread pool get initialized when the application
        ///  gets started, e.g. in `main()`.  Threads will be spawned and
        ///  initialized.  They will immediately wait on a condition variable
        ///  until a task is ready.
        ///
        Scheduler(int numThreads);

        void
        pause();

        void
        resume();

        ///
        ///  For illustration purposes it is nice to identify task by a string.
        ///  In our example these strings happen to be Latex code.  So we are
        ///  able to dump the scheduled directed graph and generate nice
        ///  pictures.
        ///
        ///  Adding a task identified by `keyStr` can depend on a list of other
        ///  tasks (`pre`) that must have been completed before it can be
        ///  started.
        ///
        void
        add(std::string keyStr, Task task, std::vector<std::string> pre= {});

        ///
        ///  In some cases it is convenient to add only some dependencies, i.e.
        ///  add an arc to the directed graph.
        ///
        ///  Note that if you want to add '`pre` must be done before `succ` gets
        ///  started' it is necessary that `succ` has not been scheduled
        ///  already.  Otherwise task `succ` could already be running.
        ///
        ///  Again tasks `pre` and `succ` are identified by string.
        ///
        void
        addArc(std::string preStr, std::string succStr);

        void
        dumpGraph(const char *filename);

        ///
        ///  Once all tasks are scheduled we wait until all tasks are completed.
        ///
        int
        join();

        ///
        ///  If a task calls `abort` then no new tasks will be stared. `join`
        ///  will wait until all running tasks are completed and returns
        ///  `errorCode`.
        ///
        void
        abort(int errorCode);

    private:

        std::pair<Task, Key>
        getTask();

        void
        runThread();

        void
        taskDone(Key key);

        Key
        getKey(std::string keyStr);

        int                                         numThreads;
        int                                         numThreadsBusy;
        int                                         errorCode;
        bool                                        stop;

        std::list<Key>                              ready;
        std::unordered_map<std::string, Key>        keyString;
        std::unordered_map<Key, std::string>        tag;
        std::unordered_map<Key, Task>               scheduled;
        std::unordered_map<Key, long>               predecessorCount;
        std::unordered_map<Key, std::list<Key>>     successor;

        std::condition_variable                     cv_scheduler;
        std::condition_variable                     cv_thread;
        std::mutex                                  mtx;
};

} // namespace flens

#endif // LU_SCHEDULER_H
