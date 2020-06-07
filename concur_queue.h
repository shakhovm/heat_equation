#ifndef MPI_LAB_CONCUR_QUEUE_H
#define MPI_LAB_CONCUR_QUEUE_H

#include <unordered_map>
#include <string>
#include <deque>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include "array2d.h"

typedef std::pair<Array2D, size_t> SEGMENT;

template<class T>
class concur_queue {
private:
    std::deque<T> queue_;
    mutable std::mutex mutex_;
    std::condition_variable cv_;
public:
    concur_queue();

    void push(T d);
    T pop();
    size_t get_size();
};


#endif //MPI_LAB_CONCUR_QUEUE_H
