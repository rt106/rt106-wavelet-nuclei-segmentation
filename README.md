# WaveletNucleiSegmentation

_Copyright (c) General Electric Company, 2017.  All rights reserved._

This is an example analytic for Rt 106. This analytic demonstrates how to use the rt106-algorithm-sdk around an [ITK](http://www.itk.org) algorithm.

### Docker container

Two Docker images are used to package and deploy the analytic.  The first is a ```dev``` Docker image that is responsible for compiling the C++ code.  The ```dev``` contains a full build environment (compilers). Building the ```dev``` image also compiles our algorithm code. After the ```dev``` image is built, a call to ```docker run``` emits a tarball containing the artifacts from the build.  Here, The second Docker image ```ops``` installs the aforementioned tarball into a runtime container.  The two container strategy allows for the runtime ```ops``` container to be small.

To build the ```dev``` container

```sh
$ cd docker
$ docker build -t rt106-dev/rt106-wavelet-nuclei-segmentation --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy=$no_proxy dev
```

To build the ```ops``` container, we first need to move the artifact from the ```dev``` container
into place to build the ```ops``` container. If we run the ```dev``` container just built, it will
emit the artifact as a tarball

```sh
$ cd docker
$ docker run rt106-dev/rt106-wavelet-nuclei-segmentation > ops/rt106-wavelet-nuclei-segmentation.tar.gz
```

To build the ```ops``` container
```sh
$ cd docker
$ docker build -t rt106/rt106-wavelet-nuclei-segmentation --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy=$no_proxy ops
```
