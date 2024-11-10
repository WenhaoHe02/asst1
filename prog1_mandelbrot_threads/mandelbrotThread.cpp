#include <stdio.h>
#include <thread>
#include <chrono>
#include <iostream>
#include "CycleTimer.h"
typedef struct
{
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int *output;
    int threadId;
    int numThreads;
} WorkerArgs;

extern void mandelbrotSerial(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int numRows,
    int maxIterations,
    int output[]);

//
// workerThreadStart --
//
// Thread entrypoint.

void workerThreadStart(WorkerArgs *const args)
{
    constexpr int stride = 2;
    for (size_t i = args->threadId * stride; i < args->height - (stride - 1); i += stride * args->numThreads)
    {
        mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height,
                         i, stride, args->maxIterations, args->output);
    }
}

// void workerThreadStart(WorkerArgs *const args)
// {
//     int half_height = 6 * args->height / 7;
//     auto start = std::chrono::high_resolution_clock::now();
//     int stride = half_height / args->numThreads;
//     int start_row = args->threadId * stride;
//     int total_row = 0;
//     if (args->threadId != args->numThreads - 2 && args->threadId != args->numThreads - 1)
//     {
//         total_row = stride;
//     }
//     else if (args->threadId == args->numThreads - 2)
//     {
//         total_row = half_height - start_row;
//     }
//     else
//     {
//         start_row = half_height;
//         total_row = args->height - start_row;
//     }

//     mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height,
//                      start_row, total_row, args->maxIterations, args->output);
//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double, std::milli> duration = end - start;
//     std::cout << "线程耗时: " << duration.count() << " 毫秒\n";
// }

// void workerThreadStart(WorkerArgs *const args)
// {
//     int thread_num = args->numThreads;
//     int stride1 = args->height / thread_num / 2;
//     int stride2 = args->height / thread_num - stride1;
//     int start_row1 = args->threadId * stride1;
//     int start_row2 = args->height - (args->threadId + 1) * stride2;
//     int total_row = 0;
//     int rest_len = args->height - (thread_num - 2) * (stride1 + stride2);

//     if (args->threadId < args->numThreads - 2) // 非最后两个线程
//     {
//         total_row = stride1;
//         mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height, start_row1, total_row, args->maxIterations, args->output);
//         total_row = stride2;
//         mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height, start_row2, total_row, args->maxIterations, args->output);
//     }
//     else if (args->threadId == args->numThreads - 2)
//     {
//         int split_row = start_row1;
//         total_row = rest_len / 2;
//         mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height, split_row, total_row, args->maxIterations, args->output);
//     }
//     else if (args->threadId == args->numThreads - 1)
//     {
//         int split_row = start_row1 - stride1 + rest_len / 2;
//         total_row = rest_len - rest_len / 2;
//         mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height, split_row, total_row, args->maxIterations, args->output);
//     }
// }

// void workerThreadStart(WorkerArgs *const args)
// {

//     // TODO FOR CS149 STUDENTS: Implement the body of the worker
//     // thread here. Each thread should make a call to mandelbrotSerial()
//     // to compute a part of the output image.  For example, in a
//     // program that uses two threads, thread 0 could compute the top
//     // half of the image and thread 1 could compute the bottom half.
//     // auto start = std::chrono::high_resolution_clock::now();
//     int thread_num = args->numThreads;
//     int stride1 = args->height / thread_num / 2;
//     int stride2 = args->height / thread_num - stride1;
//     int start_row1 = args->threadId * stride1;
//     int start_row2 = args->height - (args->threadId + 1) * stride2;
//     int total_row = 0;
//     if (args->threadId != args->numThreads - 1)
//     {
//         total_row = stride1;
//         mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height, start_row1, total_row, args->maxIterations, args->output);
//         total_row = stride2;
//         mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height, start_row2, total_row, args->maxIterations, args->output);
//     }
//     else
//     {
//         total_row = start_row2 + stride2 - start_row1;
//         mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height, start_row1, total_row, args->maxIterations, args->output);

//     }
//     // auto end = std::chrono::high_resolution_clock::now();
//     // std::chrono::duration<double, std::milli> duration = end - start;
//     // std::cout << "线程耗时: " << duration.count() << " 毫秒\n";
// }

//
// MandelbrotThread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Threads of execution are created by spawning std::threads.
void mandelbrotThread(
    int numThreads,
    float x0, float y0, float x1, float y1,
    int width, int height,
    int maxIterations, int output[])
{
    static constexpr int MAX_THREADS = 16;

    if (numThreads > MAX_THREADS)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
        exit(1);
    }

    // Creates thread objects that do not yet represent a thread.
    std::thread workers[MAX_THREADS];
    WorkerArgs args[MAX_THREADS];

    for (int i = 0; i < numThreads; i++)
    {

        // TODO FOR CS149 STUDENTS: You may or may not wish to modify
        // the per-thread arguments here.  The code below copies the
        // same arguments for each thread
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = maxIterations;
        args[i].numThreads = numThreads;
        args[i].output = output;

        args[i].threadId = i;
    }

    // Spawn the worker threads.  Note that only numThreads-1 std::threads
    // are created and the main application thread is used as a worker
    // as well.
    for (int i = 1; i < numThreads; i++)
    {
        workers[i] = std::thread(workerThreadStart, &args[i]);
    }

    workerThreadStart(&args[0]);

    // join worker threads
    for (int i = 1; i < numThreads; i++)
    {
        workers[i].join();
    }
}
