#include "EDCL_ResourceManager.h"
#include "EDCL.h"
#include "polybench.h"
#include "Stg_Sequential.h"
#include "EDCL_Scheduler.h"

#define BENCH_KERNEL syr

static EDCL *edcl;
static EDCL_Scheduler *scheduler;
static SequentialImproved *seq;
static float lambda0;

float getVoltageV() {
  return std::stof(getFileContentsStr("/sys/class/power_supply/battery/voltage_now")) / 1000000.0f;
}

float getDischargeAh() {
  return std::stof(exec("dumpsys batterystats | grep Discharge: | cut -d':' -f2 | cut -d' ' -f2")) * 0.001f;
}

#define RUN_TIME 240

float calLambda(EXECUTE_DEVICE device, EDKernel kernel) {
  scheduler->clearKernels();
  seq->setExecutionDevice(device);
  scheduler->submitKernels(kernel);
  scheduler->executeKernels(*seq);
  int loopNum = ceil(RUN_TIME / seq->executionTime);
  scheduler->clearKernels();
  for (int kI = 0; kI < loopNum; ++kI) {
	scheduler->submitKernels(edcl->createKernelCopy(kernel));
  }
  float U_before = getVoltageV();
  std::vector<float> current;
  std::vector<float> currentI0;
  float DV;
  debug_printf("Recharging...", "\n");
  sleep(RUN_TIME);
  debug_printf("Running...", "\n");
  { /* NO CHARGING ZONE */
	exec("disableUSBCharge.sh");
	// get I0
	scheduler->getResourceManager()->startCurrentRecording(currentI0);
	sleep(60);
	scheduler->getResourceManager()->stopCurrentRecording();
	scheduler->getResourceManager()->startCurrentRecording(current);
	DV = getDischargeAh();
	scheduler->executeKernels(*seq);
	DV = getDischargeAh() - DV;
	exec("enableUSBCharge.sh");
  }
  scheduler->getResourceManager()->stopCurrentRecording();
  float i0 = 0;
  for (auto &i : currentI0) {
	i0 += i;
  }
  i0 /= currentI0.size();

  float I = 0;
  for (auto &i : current) {
	I += i;
  }
  I /= current.size();
  float U_after = getVoltageV();
  auto t = (float)seq->executionTime;
  debug_printf(__func__, "DV: %f, U_before: %f, U_after: %f, I: %f, I0: %f, t: %f\n", DV, U_before, U_after, I, i0, t);
  float U = (U_before + U_after) / 2;
  float lambda = 0;
  float f0 = scheduler->getResourceManager()->getCPUCurrentFrequency(0) / 1000000.0f;
  float f1 = scheduler->getResourceManager()->getCPUCurrentFrequency(7) / 1000000.0f;
  float f2 = scheduler->getResourceManager()->getGPUCurrentFrequency() / 1000.0f;
  if (device == L_SD1 || device == L_SD2 || device == L_SD4) {
	if (device == L_SD1) {
	  lambda0 = (DV * U) / ((I * I - i0 * i0) * f0 * t);
	  lambda = lambda0;
	  debug_printf(EXECUTE_DEVICE_ToString(device), "lambda: %f\n", lambda0);
	  flushed_printf("csv_lambda0:%f, %f, %f, %f, %f, %f\n", DV, U, I, i0, f0, t);
	} else {
	  lambda = (DV * U + i0 * i0 * lambda0 * f0 * t) / (I * I * f0 * t);
	  debug_printf(EXECUTE_DEVICE_ToString(device), "lambda: %f\n", lambda);
	  flushed_printf("csv_lambda%s:%f, %f, %f, %f, %f, %f, %f, %f\n",
					 EXECUTE_DEVICE_ToString(device), DV, U, i0, lambda0, f0, t, I, f0);
	}
  } else if (device == B_SD1 || device == B_SD2 || device == B_SD4) {
	lambda = (DV * U + i0 * i0 * lambda0 * f0 * t) / (I * I * f1 * t);
	debug_printf(EXECUTE_DEVICE_ToString(device), "lambda: %f\n", lambda);
	flushed_printf("csv_lambda%s:%f, %f, %f, %f, %f, %f, %f, %f\n",
				   EXECUTE_DEVICE_ToString(device), DV, U, i0, lambda0, f0, t, I, f1);
  } else if (device == GPU) {
	lambda = (DV * U + i0 * i0 * lambda0 * f0 * t) / (I * I * f2 * t);
	debug_printf(EXECUTE_DEVICE_ToString(device), "lambda: %f\n", lambda);
	flushed_printf("csv_lambda%s:%f, %f, %f, %f, %f, %f, %f, %f\n",
				   EXECUTE_DEVICE_ToString(device), DV, U, i0, lambda0, f0, t, I, f2);
  }
  flushed_printf("Result %s lambda: %f\n", EXECUTE_DEVICE_ToString(device), lambda);
  return lambda;
}

int main() {
  edcl = new EDCL();
  scheduler = EDCL_Scheduler::getInstance(edcl);
  seq = new SequentialImproved(edcl);
  debug_printf("Voltage:", "%f v\n", getVoltageV());
  debug_printf("Discharge value:", "%f Ah\n", getDischargeAh());
  BENCH_KERNEL bench(edcl);
  auto aks = bench.createBenchAKS();
  auto rm = scheduler->getResourceManager();

  for (int loop = 0; loop < 5; ++loop) {
	for (int freq = 0; freq < 4; ++freq) {
	  flushed_printf("Frequency index%d\n", freq);
	  rm->setCPUFrequency(0, freq);
	  rm->setCPUFrequency(6, freq);
	  rm->setGPUFrequency(freq);
	  calLambda(L_SD1, aks->kernels[0]);
	  calLambda(L_SD2, aks->kernels[0]);
	  calLambda(L_SD4, aks->kernels[0]);
	  calLambda(B_SD1, aks->kernels[0]);
	  calLambda(B_SD2, aks->kernels[0]);
	  calLambda(B_SD4, aks->kernels[0]);
	  calLambda(GPU, aks->kernels[0]);
	}
  }

  rm->setSoCMaxFrequency();

  delete seq;
  delete edcl;
  return 0;
}