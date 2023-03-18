/*
**  CRC.H - header file for SNIPPETS CRC and checksum functions
*/

#ifndef CRC__H
#define CRC__H

#include <stdlib.h>           /* For size_t                           */
#include <stdint.h>           /* For uint8_t, uint16_t, uint32_t      */
#include <stdbool.h>          /* For bool, true, false                */

/*
**  File: ARCCRC16.C
*/

void init_crc_table(void);
uint16_t crc_calc(uint16_t crc, char *buf, unsigned nbytes);
void do_file(char *fn);

/*
**  File: CRC-16.C
*/

uint16_t crc16(char *data_p, uint16_t length);

/*
**  File: CRC-16F.C
*/

uint16_t updcrc(uint16_t icrc, uint8_t *icp, size_t icnt);

/*
**  File: CRC_32.C
*/

#define UPDC32(octet,crc) (crc_32_tab[((crc)\
     ^ ((uint8_t)octet)) & 0xff] ^ ((crc) >> 8))

uint32_t updateCRC32(unsigned char ch, uint32_t crc);
bool crc32file(char *name, uint32_t *crc, long *charcnt);
uint32_t crc32buf(char *buf, size_t len);

/*
**  File: CHECKSUM.C
*/

unsigned checksum(void *buffer, size_t len, unsigned int seed);

/*
**  File: CHECKEXE.C
*/

void checkexe(char *fname);



#endif /* CRC__H */
