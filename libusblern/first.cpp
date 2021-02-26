#include <errno.h>
#include <signal.h> 
#include <string.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <getopt.h>
#include <unistd.h>
#include <libusb-1.0/libusb.h> 
#include <iostream>

using namespace std;

#define VENDOR_ID 0xffff
#define PRODUCT_ID 0x0035

const static int PACKET_INT_LEN=2;
const static int INTERFACE=0; 
const static int ENDPOINT_INT_IN=0x00; /* endpoint 0x81 address for IN */ 
const static int ENDPOINT_INT_OUT=0x00; /* endpoint 0 address for OUT */ 
const static int TIMEOUT=5000; /* timeout in ms */ 

void brute_force(struct libusb_device_handle* devh){
  
}


int main(int argc,char** argv)
{ 
    int r;
    struct libusb_device_handle *devh = NULL; 
    unsigned char *data = new unsigned char[4]; //data to write
    int test = 0;

    /* Init USB */
    r = libusb_init(NULL); 
    if (r < 0) 
        fprintf(stderr, "Failed to initialise libusb\n"),exit(1); 
  
    devh = libusb_open_device_with_vid_pid(NULL, VENDOR_ID, PRODUCT_ID);

    if (!devh) { 
        fprintf(stderr, "Could not find/open LVR Generic HID device\n"); 
        goto out; 
    } 
    fprintf(stdout, "Successfully found the RFID R/W device\n"); 

    libusb_detach_kernel_driver(devh, 0);
    r = libusb_set_configuration(devh, 1); 
    if (r < 0) { 
        fprintf(stderr, "libusb_set_configuration error %d\n", r); 
        goto out; 
    } 

    r = libusb_claim_interface(devh, 0); 
    if (r < 0) { 
        fprintf(stderr, "libusb_claim_interface error %d\n", r); 
        goto out; 
    }
    data[0]='a';data[1]='b';data[2]='c';data[3]='d'; //some dummy values
    int actual; //used to find out how many bytes were written

    r = libusb_claim_interface(devh, 0); //claim interface 0 (the first) of device (mine had jsut 1)
    if(r < 0) {
        cout<<"Cannot Claim Interface"<<endl;
        return 1;
    }
    while(test < 128){
        r = libusb_bulk_transfer(devh, (test| LIBUSB_ENDPOINT_OUT), data, 4, &actual, 0); 
        if(!r) printf("%d %d\n", r, actual);
        test++;
    }
    cout<<"Claimed Interface"<<endl;

/*     libusb_device *dev;
    dev = libusb_get_device(devh);
    struct libusb_config_descriptor *config;
    r = libusb_get_config_descriptor(dev, 0, &config);
    if (r < 0) {
        fprintf(stderr, "libusb_get_config_descriptor error %d\n", r);
        goto out;
    }

    unsigned char tmp[1];
    r = libusb_control_transfer(devh,LIBUSB_ENDPOINT_IN | LIBUSB_REQUEST_TYPE_STANDARD | LIBUSB_RECIPIENT_INTERFACE,LIBUSB_REQUEST_GET_DESCRIPTOR,(LIBUSB_DT_REPORT << 8),0,tmp,1,TIMEOUT);
    libusb_reset_device(devh);
    sleep(2);  */

    brute_force(devh);
 
    libusb_release_interface(devh, 0); 
out: 
    libusb_close(devh); 
    libusb_exit(NULL); 
exit:
    return r >= 0 ? r : -r; 
} 