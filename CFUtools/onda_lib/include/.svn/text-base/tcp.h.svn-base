/*
** tcp_lib.h
** 
** Made by Morten Fischer Rasmussen
** Login   <mf@mf-black>
** 
** Started on  Wed Jul 30 20:21:19 2008 Morten Fischer Rasmussen
** Last update Wed Jul 30 20:21:19 2008 Morten Fischer Rasmussen
*/

#ifndef   	TCP_H_
# define   	TCP_H_

#define ONDA_NOT_CONNECTED 0
#define ONDA_CONNECTED     1

/* struct contraining all connection information */
typedef struct connection 
{
  char* name;
  int buf_size;
  char *rx_buf;
  char *tx_buf;
  char* addr;
  unsigned int port;
  int connected;
  int sock;
}  connection_t ;

/* typedef struct connection Connection; */


/* prototypes */
extern int tcp_init (const char* host_addr, const int host_port);
extern int tcp_term (void);
extern void tcp_auto_term (void);
extern int tcp_tx (const char* send_str);
extern int tcp_rx (char** data_ptr);
extern int tcp_query (const char* send_str, char** data_ptr);
//extern int tcp_clear_rx_buff (void);



/* error numbers */

//const char *ERR_RCV_TIMEOUT_STR ="Connection timed out. No data was received.\n";
enum {ERR_RCV_TIMEOUT  = -10};
enum {ERR_RCV_SYS_ERR = -11};
enum {ERR_RCV_NO_DATA = -12};
enum {ERR_RCV_TIMOUT2  = -13};


#endif 	    /* !TCP_H_ */
