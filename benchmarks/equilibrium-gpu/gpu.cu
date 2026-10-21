#include <cuda_runtime.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <vector>
#include "generated/ideal_map_gpu.cuh"

static void checked(cudaError_t rc) {
    if (rc != cudaSuccess) throw std::runtime_error(cudaGetErrorString(rc));
}
using Clock = std::chrono::steady_clock;
static double seconds(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now()-t).count();
}

template<int precision, int action, typename Input, typename Output>
__global__ void evaluate(int n, const Input* x, Output* y) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;
    if constexpr (precision == 0 && action == 0) {
        ideal_hessian_f64(
            x[i + 0*n], x[i + 1*n], x[i + 2*n], x[i + 3*n],
            x[i + 4*n], x[i + 5*n], x[i + 6*n], x[i + 7*n],
            x[i + 8*n], x[i + 9*n], x[i + 10*n], x[i + 11*n],
            x[i + 12*n], x[i + 13*n], x[i + 14*n], &y[i + 0*n],
            &y[i + 1*n], &y[i + 2*n], &y[i + 3*n], &y[i + 4*n],
            &y[i + 5*n], &y[i + 6*n], &y[i + 7*n], &y[i + 8*n],
            &y[i + 9*n], &y[i + 10*n], &y[i + 11*n], &y[i + 12*n],
            &y[i + 13*n], &y[i + 14*n], &y[i + 15*n], &y[i + 16*n],
            &y[i + 17*n], &y[i + 18*n], &y[i + 19*n], &y[i + 20*n],
            &y[i + 21*n], &y[i + 22*n], &y[i + 23*n], &y[i + 24*n],
            &y[i + 25*n], &y[i + 26*n], &y[i + 27*n], &y[i + 28*n],
            &y[i + 29*n], &y[i + 30*n], &y[i + 31*n], &y[i + 32*n],
            &y[i + 33*n], &y[i + 34*n], &y[i + 35*n], &y[i + 36*n],
            &y[i + 37*n], &y[i + 38*n], &y[i + 39*n], &y[i + 40*n],
            &y[i + 41*n], &y[i + 42*n], &y[i + 43*n], &y[i + 44*n]);
    }
    if constexpr (precision == 1 && action == 0) {
        ideal_hessian_f32(
            x[i + 0*n], x[i + 1*n], x[i + 2*n], x[i + 3*n],
            x[i + 4*n], x[i + 5*n], x[i + 6*n], x[i + 7*n],
            x[i + 8*n], x[i + 9*n], x[i + 10*n], x[i + 11*n],
            x[i + 12*n], x[i + 13*n], x[i + 14*n], &y[i + 0*n],
            &y[i + 1*n], &y[i + 2*n], &y[i + 3*n], &y[i + 4*n],
            &y[i + 5*n], &y[i + 6*n], &y[i + 7*n], &y[i + 8*n],
            &y[i + 9*n], &y[i + 10*n], &y[i + 11*n], &y[i + 12*n],
            &y[i + 13*n], &y[i + 14*n], &y[i + 15*n], &y[i + 16*n],
            &y[i + 17*n], &y[i + 18*n], &y[i + 19*n], &y[i + 20*n],
            &y[i + 21*n], &y[i + 22*n], &y[i + 23*n], &y[i + 24*n],
            &y[i + 25*n], &y[i + 26*n], &y[i + 27*n], &y[i + 28*n],
            &y[i + 29*n], &y[i + 30*n], &y[i + 31*n], &y[i + 32*n],
            &y[i + 33*n], &y[i + 34*n], &y[i + 35*n], &y[i + 36*n],
            &y[i + 37*n], &y[i + 38*n], &y[i + 39*n], &y[i + 40*n],
            &y[i + 41*n], &y[i + 42*n], &y[i + 43*n], &y[i + 44*n]);
    }
    if constexpr (precision == 2 && action == 0) {
        ideal_hessian_mixed(
            x[i + 0*n], x[i + 1*n], x[i + 2*n], x[i + 3*n],
            x[i + 4*n], x[i + 5*n], x[i + 6*n], x[i + 7*n],
            x[i + 8*n], x[i + 9*n], x[i + 10*n], x[i + 11*n],
            x[i + 12*n], x[i + 13*n], x[i + 14*n], &y[i + 0*n],
            &y[i + 1*n], &y[i + 2*n], &y[i + 3*n], &y[i + 4*n],
            &y[i + 5*n], &y[i + 6*n], &y[i + 7*n], &y[i + 8*n],
            &y[i + 9*n], &y[i + 10*n], &y[i + 11*n], &y[i + 12*n],
            &y[i + 13*n], &y[i + 14*n], &y[i + 15*n], &y[i + 16*n],
            &y[i + 17*n], &y[i + 18*n], &y[i + 19*n], &y[i + 20*n],
            &y[i + 21*n], &y[i + 22*n], &y[i + 23*n], &y[i + 24*n],
            &y[i + 25*n], &y[i + 26*n], &y[i + 27*n], &y[i + 28*n],
            &y[i + 29*n], &y[i + 30*n], &y[i + 31*n], &y[i + 32*n],
            &y[i + 33*n], &y[i + 34*n], &y[i + 35*n], &y[i + 36*n],
            &y[i + 37*n], &y[i + 38*n], &y[i + 39*n], &y[i + 40*n],
            &y[i + 41*n], &y[i + 42*n], &y[i + 43*n], &y[i + 44*n]);
    }
    if constexpr (precision == 0 && action == 1) {
        ideal_hvp_f64(
            x[i + 0*n], x[i + 1*n], x[i + 2*n], x[i + 3*n],
            x[i + 4*n], x[i + 5*n], x[i + 6*n], x[i + 7*n],
            x[i + 8*n], x[i + 9*n], x[i + 10*n], x[i + 11*n],
            x[i + 12*n], x[i + 13*n], x[i + 14*n], x[i + 15*n],
            x[i + 16*n], x[i + 17*n], x[i + 18*n], x[i + 19*n],
            x[i + 20*n], x[i + 21*n], x[i + 22*n], x[i + 23*n],
            &y[i + 0*n], &y[i + 1*n], &y[i + 2*n], &y[i + 3*n],
            &y[i + 4*n], &y[i + 5*n], &y[i + 6*n], &y[i + 7*n],
            &y[i + 8*n]);
    }
    if constexpr (precision == 1 && action == 1) {
        ideal_hvp_f32(
            x[i + 0*n], x[i + 1*n], x[i + 2*n], x[i + 3*n],
            x[i + 4*n], x[i + 5*n], x[i + 6*n], x[i + 7*n],
            x[i + 8*n], x[i + 9*n], x[i + 10*n], x[i + 11*n],
            x[i + 12*n], x[i + 13*n], x[i + 14*n], x[i + 15*n],
            x[i + 16*n], x[i + 17*n], x[i + 18*n], x[i + 19*n],
            x[i + 20*n], x[i + 21*n], x[i + 22*n], x[i + 23*n],
            &y[i + 0*n], &y[i + 1*n], &y[i + 2*n], &y[i + 3*n],
            &y[i + 4*n], &y[i + 5*n], &y[i + 6*n], &y[i + 7*n],
            &y[i + 8*n]);
    }
    if constexpr (precision == 2 && action == 1) {
        ideal_hvp_mixed(
            x[i + 0*n], x[i + 1*n], x[i + 2*n], x[i + 3*n],
            x[i + 4*n], x[i + 5*n], x[i + 6*n], x[i + 7*n],
            x[i + 8*n], x[i + 9*n], x[i + 10*n], x[i + 11*n],
            x[i + 12*n], x[i + 13*n], x[i + 14*n], x[i + 15*n],
            x[i + 16*n], x[i + 17*n], x[i + 18*n], x[i + 19*n],
            x[i + 20*n], x[i + 21*n], x[i + 22*n], x[i + 23*n],
            &y[i + 0*n], &y[i + 1*n], &y[i + 2*n], &y[i + 3*n],
            &y[i + 4*n], &y[i + 5*n], &y[i + 6*n], &y[i + 7*n],
            &y[i + 8*n]);
    }
}

struct Workspace {
    void *input=nullptr, *output=nullptr;
    cudaEvent_t start=nullptr, stop=nullptr;
    ~Workspace() {
        if (input) cudaFree(input);
        if (output) cudaFree(output);
        if (start) cudaEventDestroy(start);
        if (stop) cudaEventDestroy(stop);
    }
};

template<int precision, int action, typename Input, typename Output>
static void run(int n, int repeats, const double* x, double* y, double* metrics) {
    const auto cold = Clock::now();
    const int nout = action == 0 ? 45 : 9;
    Workspace w;
    const size_t ib = size_t(n)*24*sizeof(Input), ob = size_t(n)*nout*sizeof(Output);
    checked(cudaMalloc(&w.input, ib));
    checked(cudaMalloc(&w.output, ob));
    checked(cudaEventCreate(&w.start));
    checked(cudaEventCreate(&w.stop));
    std::vector<Input> packed(size_t(n)*24);
    std::copy(x, x+packed.size(), packed.begin());
    auto t = Clock::now();
    checked(cudaMemcpy(w.input, packed.data(), ib, cudaMemcpyHostToDevice));
    metrics[1] = seconds(t);
    auto launch = [&]() {
        evaluate<precision,action><<<(n+127)/128,128>>>(n,
            static_cast<Input*>(w.input), static_cast<Output*>(w.output));
        checked(cudaGetLastError());
    };
    std::vector<Output> output(size_t(n)*nout);
    launch();
    checked(cudaMemcpy(output.data(),w.output,ob,cudaMemcpyDeviceToHost));
    metrics[0]=seconds(cold);
    std::vector<double> samples(repeats);
    for (int r=0; r<repeats; ++r) {
        checked(cudaEventRecord(w.start));
        launch();
        checked(cudaEventRecord(w.stop));
        checked(cudaEventSynchronize(w.stop));
        float ms;
        checked(cudaEventElapsedTime(&ms,w.start,w.stop));
        samples[r]=ms*1e-3;
    }
    std::sort(samples.begin(), samples.end());
    metrics[2]=samples.front(); metrics[3]=samples[repeats/2]; metrics[4]=samples.back();
    t = Clock::now();
    checked(cudaMemcpy(output.data(),w.output,ob,cudaMemcpyDeviceToHost));
    metrics[5]=seconds(t);
    std::copy(output.begin(),output.end(),y);
    metrics[6]=double(ib+ob);
}

extern "C" int kin6d_gpu_benchmark(int n, int repeats, int precision, int action,
                                   const double* x, double* y, double* metrics) {
    try {
        if (n<1 || repeats<1) return 1;
        if (action == 0) {
            if (precision == 0) run<0,0,double,double>(n,repeats,x,y,metrics);
            else if (precision == 1) run<1,0,float,float>(n,repeats,x,y,metrics);
            else if (precision == 2) run<2,0,float,double>(n,repeats,x,y,metrics);
            else return 1;
        } else if (action == 1) {
            if (precision == 0) run<0,1,double,double>(n,repeats,x,y,metrics);
            else if (precision == 1) run<1,1,float,float>(n,repeats,x,y,metrics);
            else if (precision == 2) run<2,1,float,double>(n,repeats,x,y,metrics);
            else return 1;
        } else return 1;
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr,"CUDA prototype: %s\n", e.what());
        return 1;
    }
}
