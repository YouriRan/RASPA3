module;

#include <cstddef>
#include <optional>
#include <string>
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
  #include <OpenCL/cl.h>
#elif _WIN32
  #include <CL/cl.h>
#else
  #include <CL/opencl.h>
#endif


export module opencl;

export namespace OpenCL
{
  extern std::optional<cl_context> clContext;
  extern std::optional<cl_device_id> clDeviceId;
  extern std::optional<cl_command_queue> clCommandQueue;

  void initialize(void);
  std::optional<cl_device_id> bestOpenCLDevice(cl_device_type device_type);
  std::string printBestOpenCLDevice();
  bool supportsImageFormatCapabilities(cl_context &trial_clContext, cl_device_id &trial_clDeviceId);
}

