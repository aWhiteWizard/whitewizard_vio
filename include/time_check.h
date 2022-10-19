#pragma once

#include <ctime>
#include <cstdlib>
#include <chrono>

class timecheck
{
  public:
    timecheck()
    {
        check_time();
    }

    void check_time()
    {
        start = std::chrono::system_clock::now();
    }

    double sum_time()
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count() * 1000;
    }

  private:
    std::chrono::time_point<std::chrono::system_clock> start, end;
};
