#ifndef DEVICECUDA_H
#define DEVICECUDA_H

#define CUDA_EVENT_ID_POLL_REMOTE 98
#define CUDA_TRACE_POLL_REMOTE \
  traceUserEvent(CUDA_EVENT_ID_POLL_REMOTE)
#define CUDA_EVENT_ID_POLL_LOCAL 99
#define CUDA_TRACE_POLL_LOCAL \
  traceUserEvent(CUDA_EVENT_ID_POLL_LOCAL)
#define CUDA_EVENT_ID_BASE 100
#define CUDA_TRACE_REMOTE(START,END) \
  do { int dev; cudaGetDevice(&dev); traceUserBracketEvent( \
       CUDA_EVENT_ID_BASE + 2 * dev, START, END); } while (0)
#define CUDA_TRACE_LOCAL(START,END) \
  do { int dev; cudaGetDevice(&dev); traceUserBracketEvent( \
       CUDA_EVENT_ID_BASE + 2 * dev + 1, START, END); } while (0)

#ifdef WIN32
#define __thread __declspec(thread)
#endif

//
// Class that handles PE <=> CUDA device mapping
//
class DeviceCUDA {

private:
	// Command line argument settings
	char *devicelist;
	int usedevicelist;
	int ignoresharing;
	int mergegrids;
	int nomergegrids;
	int nostreaming;

	// True when GPU is shared between PEs
	bool sharedGpu;
	// Index of next GPU sharing this GPU
	int nextPeSharingGpu;
	// Index of the master PE for this GPU
	int masterPe;
	// Number of PEs that share this GPU
	int numPesSharingDevice;
	// List of PEs that share this GPU
	int *pesSharingDevice;
	// True when what???
	int gpuIsMine;


	void register_user_events();

public:
	DeviceCUDA();
	~DeviceCUDA();
	
	void initialize();

	bool device_shared_with_pe(int pe);
	bool one_device_per_node();

	int getNoStreaming() {return nostreaming;}
	int getNoMergeGrids() {return nomergegrids;}
	int getMergeGrids() {return mergegrids;}
	void setMergeGrids(const int val) {mergegrids = val;}

	bool getSharedGpu() {return sharedGpu;}
	int getNextPeSharingGpu() {return nextPeSharingGpu;}
	int getMasterPe() {return masterPe;}
	int getNumPesSharingDevice() {return numPesSharingDevice;}
	int getPesSharingDevice(const int i) {return pesSharingDevice[i];}

	int getGpuIsMine() {return gpuIsMine;}
	void setGpuIsMine(const int val) {gpuIsMine = val;}
};

#endif // DEVICECUDA_H
