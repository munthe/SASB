#include "mex.h"

#if defined __linux__
#include <netdb.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#define SOCKET int
#elif defined __APPLE__
#include <netdb.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#define SOCKET int
#elif defined _MSC_VER
#include <winsock2.h>
#include <ws2tcpip.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

int init_network(char* hostname, char* port_str);
void send_network(char *command);
void recv_network(char* &response, size_t buf_size);
void end_network(void);
void clear_rx_buffer(void);

//const size_t max_num_sockets=10;
SOCKET SocketFD;
bool is_socket_inited;

const size_t buffer_size=4096;

/* Input arguments besides socket index and command:
   Init: IP address, port
   Send: String containing command
   Read: None
   End: None */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t i;
  char hostname[65];
  char port_str[6];
  char command[5];
  char remote_command[buffer_size];
  char* remote_answer;
  remote_answer=(char*) calloc(buffer_size, sizeof(char));
  int ret;

  if (nrhs < 1) {
    mexErrMsgTxt("ERROR *** No arguments given!");
  }

  if (1 == mxGetString(prhs[0], command, 5)) {
    mexErrMsgTxt("ERROR *** Conversion of input to command failed");
  }
  command[4] = 0;


  /* INIT */
  if (0 == strcmp(command, "init")) {
    if (nrhs != 3) {
      mexErrMsgTxt("ERROR *** Exactly 2 arguments needed for init");
    }
    if (1 == mxGetString(prhs[1], hostname, 65)) {
      mexErrMsgTxt("ERROR *** Conversion of input to hostname failed");
    }
    if (1 == mxGetString(prhs[2], port_str, 6)) {
      mexErrMsgTxt("ERROR *** Conversion of input to port string failed");
    }
    hostname[64] = 0;
    port_str[5] = 0;
    
    ret = init_network(hostname, port_str);
    if (ret) {
      mexErrMsgTxt("ERROR *** Error initializing network connection");
    }



    /* SEND */
  } else if (0 == strcmp(command, "send")) {
    if (nrhs != 2) {
      mexErrMsgTxt("ERROR *** Exactly 1 argument needed for send");
    }
    if (!is_socket_inited) {
      mexErrMsgTxt("ERROR *** Attempting to use uninitialized socket");
    }
    if (1 == mxGetString(prhs[1], remote_command, buffer_size)) {
      mexErrMsgTxt("ERROR *** Error getting command to execute");
    }
    remote_command[buffer_size-1]=0;
    send_network(remote_command);



    /* RCV */
  } else if (0 == strcmp(command, "recv")) {
    if (nlhs != 1) {
      mexErrMsgTxt("ERROR *** Exactly 1 output produced for recv");
    }
    if (!is_socket_inited) {
      mexErrMsgTxt("ERROR *** Attempting to use uninitialized socket");
    }
    recv_network(remote_answer, buffer_size);
    remote_answer[buffer_size-1]=0;
    plhs[0]=mxCreateString(remote_answer);



    /* END */
  } else if (0 == strcmp(command, "end")) {
    if (!is_socket_inited) {
      mexErrMsgTxt("ERROR *** Attempting to use uninitialized socket");
    }
    end_network();
  } else {
    mexErrMsgTxt("ERROR *** Unknown command specified");
  }
}




int init_network(char* hostname, char* port_str)
{
  unsigned int to_send;
  unsigned char buf[2];
  unsigned int status;
  
#if defined __linux__
#elif defined WIN32
  WSADATA wsaData;
  SOCKET ListenSocket = INVALID_SOCKET;
  int iResult;

  // Initialize Winsock
  iResult = WSAStartup(MAKEWORD(2,2), &wsaData);
  if (iResult != 0) {
    mexPrintf("WSAStartup failed with error: %d\n", iResult);
    return -1;
  }
#endif

  struct addrinfo* _addrinfo;
  struct sockaddr stTCPSockAddr;

  SocketFD = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);

  if (-1 == SocketFD) {
    mexErrMsgTxt("cannot create socket");
  }
 
  memset(&stTCPSockAddr, 0, sizeof(sockaddr_in));
  if (0 != getaddrinfo(hostname, port_str, NULL, &_addrinfo)) {
    return -1;
  }
  stTCPSockAddr = *(_addrinfo->ai_addr);

  if (-1 == connect(SocketFD, &stTCPSockAddr, sizeof(sockaddr))) {
    mexErrMsgTxt("connect failed");
  }
  is_socket_inited=true;



  /* automatic close connection when MATLAB clears or exits */
  mexAtExit(end_network);

  return 0;
}




void end_network(void)
{
  if (shutdown(SocketFD, SHUT_RDWR)) {
    mexErrMsgTxt("ERROR *** Shutting down socket failed");
  }
  if (close(SocketFD)) {
    mexErrMsgTxt("ERROR *** Closing socket failed");
  }
  is_socket_inited=false;
}



void send_network(char* command)
{
  bool done = false;
  size_t n_bytes;
  size_t n_sent;
  size_t to_send=strlen(command);
  

  clear_rx_buffer(); // make sure the receive buffer is empty

  // Send command
  n_sent = 0;
  while (n_sent != to_send) {
    n_bytes = send(SocketFD, command+n_sent, to_send-n_sent, 0);
    if (n_bytes == -1) {
#if defined __linux__
      mexPrintf(strerror(errno));
#elif defined _MSC_VER
      LPVOID lpMsgBuf;
      int err_code = WSAGetLastError();
      if (FormatMessage(
              FORMAT_MESSAGE_ALLOCATE_BUFFER |
              FORMAT_MESSAGE_FROM_SYSTEM |
              FORMAT_MESSAGE_IGNORE_INSERTS,
              NULL, err_code,
              MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
              // Default language
              (LPTSTR)&lpMsgBuf, 0, NULL)) {
        mexPrintf("Error %d: %s\n", err_code, lpMsgBuf);
        LocalFree(lpMsgBuf);
      } else {
        mexPrintf("Error %d\n", err_code);
      }
#endif
      mexErrMsgTxt("ERROR *** Error sending data");
    }
    n_sent += n_bytes;
  }
}



#define RECEIVE(status, n_bytes_rcv, socket, var_ptr, length, flags)	\
  n_bytes_rcv = 0;							\
  while ((size_t) n_bytes_rcv < length) {				\
    status = recv(socket, var_ptr, length-(size_t) n_bytes_rcv, flags);	\
    if (status < 0) {							\
      mexErrMsgTxt("ERROR *** Error receiving data");			\
    }									\
    if (status == 0) {							\
      break;								\
    }									\
    n_bytes_rcv += status;						\
    if (var_ptr[n_bytes_rcv-1]=='\n') { \
      var_ptr[n_bytes_rcv-1]=0; \
      if (var_ptr[n_bytes_rcv-2]=='\r') { \
        var_ptr[n_bytes_rcv-2]=0; \
      } else { \
        mexErrMsgTxt("ERROR *** Unexpected line end"); \
      } \
      break; \
    } \
  }



void recv_network(char* &response, size_t buf_size)
{
  size_t recv_status;
  size_t n_bytes_rcv;

  RECEIVE(recv_status, n_bytes_rcv, SocketFD, response, buf_size, 0);
}



void clear_rx_buffer(void)
{
  fd_set read_fd;
  struct timeval timeout;
  int nfds = SocketFD+1;
  size_t buf_size = 256;
  char rx_buf[buf_size];
  int retval;

  /* clear the struct */
  FD_ZERO (&read_fd);
  /* add the socket to the struct/list */
  FD_SET (SocketFD, &read_fd);
  /* set timeout to 0.001 sec */
  timeout.tv_sec = 0;
  timeout.tv_usec = 1000;
  /* Wait until data is available or timeout has passed */
  retval = select (nfds, &read_fd, NULL, NULL, &timeout);

  if (retval > 0) /* data present in buffer */
    {
      /* clear buffer */
      /* data available -> receive it */
      recv (SocketFD, rx_buf, buf_size -1, 0);
    }
}
