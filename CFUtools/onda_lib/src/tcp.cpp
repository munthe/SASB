/*
** tcp.cpp
**
** This file implements the TCP connection to the Onda system.
** It was originally written in C for the SARUS project.
**
** $Id: tcp.cpp 6 2012-04-16 15:40:01Z mf $
**
** Made by (Morten Fischer Rasmussen)
** Login   <mf@mf-black>
**
** Started on  Wed Jul 30 20:22:23 2008 Morten Fischer Rasmussen
** Modified for Onda on Wed Feb 22 12:27:10 2012 Morten Fischer Rasmussen
*/

#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>  /* provides close () */
#include <netinet/in.h>
/* #include <fcntl.h> */
#include <errno.h>
#include <sys/select.h>
#include <sys/time.h>

#include "tcp.h"
#include "err_printf.h"


#define ONDA_BUFF_SIZE 9096
#define ONDA_PORT 49999


/* input/output buffers */
static char rx_buf[ONDA_BUFF_SIZE];
static char tx_buf[ONDA_BUFF_SIZE];
static char onda_addr[30];


static connection_t onda_st = {(char*)"Onda",
			       ONDA_BUFF_SIZE,
			       rx_buf,
			       tx_buf,
			       onda_addr,
			       ONDA_PORT,
			       ONDA_NOT_CONNECTED,
			       0};



/**
 *  Clears the local buffer and system buffer
 */
static int
tcp_clear_rx_buff (void)
{
    int retval;
    int nfds = onda_st.sock + 1;
    fd_set read_fd;
    struct timeval timeout;
    /* make sure pointer does not contain random value */
    // *data_ptr = NULL;

    if (onda_st.connected != ONDA_CONNECTED){
	err_printf ("Error on clearing buffer: Cannot receive on non-existing connection.\n");
	return -1;
    }

    /* clear the struct */
    FD_ZERO (&read_fd);
    /* add the socket to the struct/list */
    FD_SET (onda_st.sock, &read_fd);
    /* set timeout to 0.001 sec */
    timeout.tv_sec = 0;
    timeout.tv_usec = 1000;
    /* Wait until data is available or timeout has passed */
    retval = select (nfds, &read_fd, NULL, NULL, &timeout);

    if (retval > 0) /* data present in buffer */
	{
	    /* clear buffer */
	    /* data available -> receive it */
	    recv (onda_st.sock, onda_st.rx_buf, onda_st.buf_size -1, 0);
	}

    /* Reset local buffer */
    rx_buf [0] = '\0';

    return 0;
}





/**
 * Creates a socket that uses TCP/IP and connects the destination adr+port
 * Input:   -struct containing all connection info
 *          -pointer to  a char buffer where error messages are printed
 *
 * Returns:  0 on success
 *          -1 on failure
 */
int
tcp_init (const char *addr_host, const int port_host)
{

    if (onda_st.connected == ONDA_CONNECTED)
	{
	    err_printf("Onda: Connection already initialised.\n");
	    return -1;
	}
    sprintf(onda_st.addr, "%s", addr_host);
    onda_st.port = port_host;
    /* clear buffers */
    rx_buf[0] = '\0';
    tx_buf[0] = '\0';


    struct sockaddr_in dst_server;

    if (onda_st.connected != ONDA_NOT_CONNECTED){
	info_printf ("Onda interface err.: Can not initialize an already excisting connection.\n",
		     onda_st.name);
	return -1;
    }


    /* input validation */
    if (onda_st.port > 65535){
	info_printf ("Failed to create socket for %s: maximum value of destination port is 65535\n", onda_st.name);
	return -1;
    }
    /* We assume that destination adresses are always valid (this is not necessary true) */
    /* Create the IP/TCP socket */
    if ((onda_st.sock = socket (PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0){
	info_printf ("Failed to create socket for %s.\nSystem returned: %s.\n", onda_st.name, strerror(errno));
	return -1;
    }

    /* Construct the server sockaddr_in structure */
    memset (&dst_server, 0, sizeof (dst_server));          /* Clear struct */
    dst_server.sin_family = AF_INET;                       /* Internet/IP */
    dst_server.sin_addr.s_addr = inet_addr (onda_st.addr);   /* IP address */
    dst_server.sin_port = htons (onda_st.port);              /* server port */

    /* Establish Connection */
    if (connect (onda_st.sock, (struct sockaddr *) &dst_server, sizeof(dst_server))  < 0){
	info_printf ("Error: Failed to connect with %s on address: %s port: %i.\nSystem returned: %s.\n",
		     onda_st.name, onda_st.addr, onda_st.port, strerror(errno));
	return -1;
    }
    else
	onda_st.connected = ONDA_CONNECTED;

    return 0;
}







int
tcp_term (void)
{
    int retval;
    if (onda_st.connected != ONDA_CONNECTED){
	info_printf ("Onda interface err.: Can not close non existing connection: %s.\n",
		     onda_st.name);
	return -1;
    }

    retval = close(onda_st.sock);
    if (retval == 0){
	onda_st.connected = ONDA_NOT_CONNECTED;
	onda_st.sock = 0;
	info_printf("Onda interface: closed connection.\n");
    }else{
	info_printf ("Onda interface err.: Could not terminate connection.\nSytem returned: %s.\n", strerror(errno));
	return retval;
    }

    return 0;
}






void
tcp_auto_term (void)
{
    if (onda_st.connected == ONDA_CONNECTED)
	tcp_term();

    return;
}








/**
 * Sends raw data using the socket handle
 * Input:   -struct containing all connection info
 *          -pointer to a char buffer where error messages are printed
 *
 * Returns:  0 on success
 *          -1 on failure
 */
int
tcp_tx (const char* send_str)
{
    int echolen;
    int retval;
    /* set the maximum number of file descriptors */
    int nfds = onda_st.sock + 1;
    fd_set write_fd;
    struct timeval timeout;

    if (onda_st.connected != ONDA_CONNECTED){
	info_printf ("Error: can not send data. Connection is non initialised.\n");
	return -1;
    }

    /* Clear the receive buffers */
    retval = tcp_clear_rx_buff();
    if (retval)
	return -1;

    /* clear the struct */
    FD_ZERO (&write_fd);
    /* add the socket to the struct/list */
    FD_SET (onda_st.sock, &write_fd);

    /* set timeout to 5 sec */
    timeout.tv_sec = 5;
    timeout.tv_usec = 0;

    /* Wait until there is room in the socket send buffer or select() times out */
    retval = select (nfds, NULL, &write_fd, NULL, &timeout);

    switch (retval)
	{
	    /* timeout */
	case 0:
	    info_printf ("Unable to send data. Socket is not ready to send data.\n");
	    return -1;
	    break;

	    /* error */
	case -1:
	    info_printf ("Unable to send data.\nSystem returned: %s\n", strerror(errno));
	    return -1;
	    break;

	    /* AOK */
	default: break;
	}


    /* Send the string to the server */
    /* echolen = strlen (onda_st.tx_buf); */
    /* retval  = send (onda_st.sock, onda_st.tx_buf, echolen, 0); /\* "0" for IP-protocol *\/ */
    echolen = strlen (send_str);
    retval  = send (onda_st.sock, send_str, echolen, 0); /* "0" for IP-protocol */

    /* test for success */
    if (retval == -1)
	{
	    info_printf ("Unable to send data\nSystem returned: %s\n", strerror(errno));
	    return -1;
	}
    else if (retval != echolen)
	{
	    info_printf ("Unable to send all data. Should have sent: %i, but only %i bytes was sent.\nSystem returned: %s\n",
			 echolen, retval, strerror(errno));
	    return -1;
	}
    return 0;
}







/**
 * Receives data using the socket handle
 * Input:   -struct containing all connection info
 *          -pointer to a char buffer where error messages are printed
 *
 * Returns:  0 on success
 *          -1 on failure
 */
int
tcp_rx (char** data_ptr)
{
    int retval;
    int nfds = onda_st.sock + 1;
    fd_set read_fd;
    struct timeval timeout;
    /* make sure pointer does not contain random value */
    // *data_ptr = NULL;

    if (onda_st.connected != ONDA_CONNECTED){
	info_printf ("Onda interface err.: Can not receive on non-existing connection.\n");
	return -1;
    }


    /* clear the struct */
    FD_ZERO (&read_fd);
    /* add the socket to the struct/list */
    FD_SET (onda_st.sock, &read_fd);
    /* set timeout to 5 sec */
    timeout.tv_sec = 5;
    timeout.tv_usec = 0;
    /* Wait until data is available or timeout has passed */
    retval = select (nfds, &read_fd, NULL, NULL, &timeout);

    /* Wait a little extra time, to assure the entire string is received. */
    timeout.tv_sec = 0;
    timeout.tv_usec = 50000; /* 0.05sec */
    select (nfds, NULL, NULL, NULL, &timeout);

    switch (retval)
	{
	    /* time out */
	case 0:
	    return ERR_RCV_TIMEOUT;
	    break;

	    /* error */
	case -1:
	    info_printf ("Error: could not receive data.\nSystem returned: %s.\n", strerror(errno));
	    return ERR_RCV_SYS_ERR;
	    break;

	    /* AOK */
	default: break;
	}

    /* data available -> receive it */
    if ((retval = recv (onda_st.sock, onda_st.rx_buf, onda_st.buf_size -1, 0)) < 0){
	rx_buf [0] = '\0';
	/* error handling */
	info_printf("Error: could not receive data.\nSystem returned: %s.\n", strerror(errno));
	/* Assure null terminated string */
	return ERR_RCV_NO_DATA;
    }
    else
	/* Assure null terminated string */
	rx_buf [retval] = '\0';


    /* set pointer to buffer */
    *data_ptr = rx_buf;

    return 0;
}






/**
 * Sends AND receives data using the socket handle
 * Input:   -c-string containing the TX string
 *          -pointer to a char buffer where error messages are printed
 *
 * Returns:  0 on success
 *          -1 on failure
 */
int
tcp_query (const char* send_str, char** receive_ptr)
{
    int retval;
    retval = tcp_tx(send_str);
    if (retval)
	return retval;

    /* receive string */
    retval = tcp_rx(receive_ptr);
    if (retval)
	return retval;

    return 0;
}
