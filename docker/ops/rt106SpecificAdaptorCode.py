# Copyright (c) General Electric Company, 2017.  All rights reserved.

# Rt 106

# Wavelet-Based Nuclei Segmentation algorithm adaptor

import os, glob, uuid, time, logging, subprocess


# function: run_algorithm()
# parameters:
#   datastore: object to be used when interacting with the Rt. 106 datastore
#   context:  A JSON structure that should contain all the inputs and parameters your algorithm needs.
def run_algorithm(datastore,context):

    logging.info('run_algorithm: %r' % context)

    # cleanup the input and output directories
    for f in glob.glob('/rt106/input/*') + glob.glob('/rt106/output/*'):
        os.remove(f)

# 1.    Code for marshalling inputs.
    if (context['slide'] == "" or context['region'] == "" or context['channel'] == ""):
        status = "ERROR_INPUT_NOT_FULLY_SPECIFIED"
        return { 'result' : {}, 'status' : status }
    input_path = datastore.get_pathology_primary_path(context['slide'], context['region'], context['channel'])
    if (type(input_path) == "int" and input_path > 200):
        status = "ERROR_INPUT_FILE_NOT_FOUND"
        return { 'result' : {}, 'status' : input_path }
    input_image = 'DAPI.tiff'
    input_file = '/rt106/input/%s' % input_image
    logging.info("WaveletNucleiSegmentation: input_image is " + input_image)
    logging.info("WaveletNucleiSegmentation: input_path is " + input_path)
    instance_status = datastore.get_instance(input_path,'/rt106/input', input_image, 'tiff16')
    if (instance_status != 200):
        status = "ERROR_INPUT_FILE_NOT_FOUND"
        return { 'result' : {}, 'status' : status }

    output_path = datastore.get_pathology_result_path(context['slide'], context['region'], context['branch'], 'NucSeg')
    logging.info("WaveletNucleiSegmentation: output_path is " + output_path)
    output_image = 'SegMask.tiff'
    output_file = '/rt106/output/%s' % output_image

    # 2.    Code for calling algorithm.
    try:
        io_args = '%s %s' % (input_file,output_file)
        algo_args = '%s %s %s %s' % (context['minSize'],context['maxSize'],context['DetectionSize'],context['sensitivity'])
        run_algorithm = '/rt106/bin/itkWaveletNucleiSegmentationTest %s %s' % (io_args,algo_args)
        logging.info('run Algorithm: %r' % run_algorithm)
        subprocess.check_call(run_algorithm,shell=True)

    except subprocess.CalledProcessError, e:
        logging.error('%d - %s' % (e.returncode, e.cmd))
        status = "EXECUTION_FINISHED_ERROR"
        result_context = {}
        return { 'result' : result_context, 'status' : status }

    # 3.    Set  status as appropriate.
    status = "EXECUTION_FINISHED_SUCCESS"

    # 4.    Store results in datastore.
    response_upload = datastore.post_instance(output_path,  '/rt106/output', output_image,  'tiff16', context['force'])

    if response_upload == 403:
        status = "EXECUTION_ERROR"

    try:
        logging.info("WaveletNucleiSegmentation: response_upload is " + str(response_upload))
    except:
        logging.info("WaveletNucleiSegmentation, error logging the status")

    logging.info("WaveletNucleiSegmentation: returned full path is %s", response_upload['path'])

    # 5.    Create JSON structure containing results.
    result_context = {
        "nucleiImage" : input_path,
        "nucleiMap" : response_upload['path']
    }

    return { 'result' : result_context, 'status' : status }
