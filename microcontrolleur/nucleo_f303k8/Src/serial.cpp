/** Includes **/
/**************/
#include "serial.h"
#include <string.h>
using namespace std;

/** Constructor **/
/*****************/
Serial::Serial(UART_HandleTypeDef* serial)
{
	m_serial_interface_ptr = serial;
}

void Serial::send_number(uint32_t value)
{
	// char buf[sizeof(uint32_t)+1];
	// snprintf(buf, sizeof buf, "%lu", (unsigned long)value);
	uint8_t buf[sizeof(uint32_t)];
	uint8_t i;
	for( i=0; i < sizeof(uint32_t); i++)
	{
		buf[i] = value >> 8*(sizeof(uint32_t)-i-1);
		send_number(buf[i]);
	}
	
}

void Serial::send_number(uint16_t value)
{
	uint8_t buf[sizeof(uint16_t)];
	uint8_t i;
	for( i=0; i < sizeof(uint16_t); i++)
	{
		buf[i] = value >> 8*(sizeof(uint16_t)-i-1);
		send_number(buf[i]);
	}

}

void Serial::send_number(uint8_t value)
{
	HAL_UART_Transmit(m_serial_interface_ptr, &value, 1, 0xFFFF);
}

void Serial::send_string(char* msg)
{
	HAL_UART_Transmit(m_serial_interface_ptr, (uint8_t*)msg, strlen(msg), 0xFFFF );
}

void Serial::send_string(uint32_t value)
{
	char buf[11];
	snprintf(buf, sizeof buf, "%lu", (unsigned long)value);
	send_string(buf);
}

void Serial::send_string(uint16_t value)
{
	char buf[6];
	snprintf(buf, sizeof buf, "%u", value);
	send_string(buf);
}

void Serial::send_string(uint8_t value)
{
	char buf[4];
	snprintf(buf, sizeof buf, "%hu", value);
	send_string(buf);
}

void Serial::send_string(int32_t value)
{
	char buf[12];
	snprintf(buf, sizeof buf, "%li", value);
	send_string(buf);
}

void Serial::send_string(int16_t value)
{
	char buf[7];
	snprintf(buf, sizeof buf, "%i", value);
	send_string(buf);
}

void Serial::send_string(int8_t value)
{
	char buf[5];
	snprintf(buf, sizeof buf, "%hi", value);
	send_string(buf);
}