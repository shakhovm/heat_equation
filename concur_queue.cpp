#include "concur_queue.h"

template class concur_queue<SEGMENT>;

template<class T>
concur_queue<T>::concur_queue() { }

template<class T>
void concur_queue<T>::push(T d) {
    {
        std::lock_guard<std::mutex> lg{mutex_};
        queue_.push_back(d);
    }
    cv_.notify_one();
}

template<class T>
T concur_queue<T>::pop() {
    std::unique_lock<std::mutex> lg{mutex_};
    cv_.wait(lg, [this](){return queue_.size() != 0;});
    T d = queue_.front();
    queue_.pop_front();
    return d;
}

template<class T>
size_t concur_queue<T>::get_size() {
    std::lock_guard<std::mutex> lg{mutex_};
    return queue_.size();
}
