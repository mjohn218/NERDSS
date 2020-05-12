// this file is copied from https://www.cnblogs.com/MakeView660/p/12531566.html, used for tracing function called
#pragma once
#ifndef TRACING_H
#define TRACING_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

// Simple structure to handle function-call tracing.

// On debug builds, always build with tracing enabled unless explicitly disabled
#if defined(_DEBUG) && !defined(TRACING_DISABLED)
#define TRACING_ENABLED
#endif

// Define a preprocessor macro to help with the tracing
#ifdef TRACING_ENABLED
#define TRACE() \
    tracing::tracer _tracer_object__##__COUNTER__ { __func__, __FILE__, __LINE__ }
#else
#define TRACE() // Nothing
#endif

#ifdef TRACING_ENABLED
namespace tracing {
class tracer {
public:
    tracer() = delete; // Disallow default construction
    tracer(tracer const&) = delete; // Disallow copy construction
    tracer(tracer&&) = delete; // Disallow move construction
    tracer& operator=(tracer const&) = delete; // Disallow copy assignment
    tracer& operator=(tracer&&) = delete; // Disallow move assignment

    tracer(std::string const& fun, std::string const& file, int const line)
        : function_name { fun }
        , file_name { file }
        , line_number { line }
    {
        std::clog << "TRACE: Entering function " << function_name << " (" << file_name << ':' << line_number << ')' << std::endl;
        startTime = std::chrono::system_clock::now();
        auto startTimeFormat = std::chrono::system_clock::to_time_t(startTime);
        char charTime[24];
        std::clog << "Entering time: ";
        if (0 < strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&startTimeFormat)))
            std::clog << charTime << std::endl;
    }

    ~tracer()
    {
        std::clog << "TRACE: Leaving function " << function_name << std::endl;
        auto endTime = std::chrono::system_clock::now();
        auto endTimeFormat = std::chrono::system_clock::to_time_t(endTime);
        char charTime[24];
        std::clog << "Leaving time: ";
        if (0 < strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&endTimeFormat)))
            std::clog << charTime << std::endl;
        std::chrono::duration<double> wallTime = endTime - startTime;
        std::clog << "\tWall Time: ";
        std::cout << wallTime.count() << " seconds" << std::endl;
    }

private:
    std::string function_name;
    std::string file_name;
    int line_number;
    std::chrono::system_clock::time_point startTime {};
};
}
#endif // TRACING_ENABLED

#endif // TRACING_H