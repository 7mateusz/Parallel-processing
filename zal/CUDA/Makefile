NVCC = nvcc
NVCCFLAGS = -O3 -arch=sm_86
TARGET = cuda

all: $(TARGET)

$(TARGET): cuda.cu
	$(NVCC) $(NVCCFLAGS) -o $(TARGET) cuda.cu

clean:
	rm -f $(TARGET) config_* result_*

.PHONY: all clean
