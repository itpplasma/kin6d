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
        if constexpr (precision == 0 && action == 0)
            ideal_hessian_f64_batch<<<(n+127)/128,128>>>(n,
                static_cast<Input*>(w.input), static_cast<Output*>(w.output));
        if constexpr (precision == 0 && action == 1)
            ideal_hvp_f64_batch<<<(n+127)/128,128>>>(n,
                static_cast<Input*>(w.input), static_cast<Output*>(w.output));
        if constexpr (precision == 1 && action == 0)
            ideal_hessian_f32_batch<<<(n+127)/128,128>>>(n,
                static_cast<Input*>(w.input), static_cast<Output*>(w.output));
        if constexpr (precision == 1 && action == 1)
            ideal_hvp_f32_batch<<<(n+127)/128,128>>>(n,
                static_cast<Input*>(w.input), static_cast<Output*>(w.output));
        if constexpr (precision == 2 && action == 0)
            ideal_hessian_mixed_batch<<<(n+127)/128,128>>>(n,
                static_cast<Input*>(w.input), static_cast<Output*>(w.output));
        if constexpr (precision == 2 && action == 1)
            ideal_hvp_mixed_batch<<<(n+127)/128,128>>>(n,
                static_cast<Input*>(w.input), static_cast<Output*>(w.output));
        checked(cudaGetLastError());
    };
    std::vector<Output> output(size_t(n)*nout);
    launch();
    checked(cudaMemcpy(output.data(),w.output,ob,cudaMemcpyDeviceToHost));
    std::copy(output.begin(),output.end(),y);
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
