#ifndef __SERIAL_H__
#define __SERIAL_H__


/** Includes **/
/**************/
#include "stm32f3xx_hal.h"
#include "stm32f3xx_hal_uart.h"

/** Defines **/
/*************/
class Serial
{
	private:
		UART_HandleTypeDef* m_serial_interface_ptr;
	public:
		Serial(UART_HandleTypeDef* serial);
		void send_number(uint32_t value);
		void send_number(uint16_t value);
		void send_number(uint8_t value);

		void send_string(char* msg);
		void send_string(uint32_t value);
		void send_string(uint16_t value);
		void send_string(uint8_t value);

};




#endif //SERIAL_H