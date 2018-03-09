#ifndef __SERIAL_H__
#define __SERIAL_H__


/** Includes **/
/**************/

#include "stm32f3xx_hal.h"
#include "stm32f3xx_hal_dma.h"
#include "stm32f3xx_hal_uart.h"

// #include "stm32f3xx_hal_tim.h"


/** Defines **/
/*************/
class Serial
{
	private:
		UART_HandleTypeDef* m_serial_interface_ptr;
	public:
		Serial(UART_HandleTypeDef* serial);
		~Serial();
		void write(uint32_t value);
		void write(uint16_t value);
		void write(uint8_t value);
		void write(uint8_t* msg,uint16_t len);

		void print(char* msg);
		void print(uint32_t value);
		void print(uint16_t value);
		void print(uint8_t value);
		void print(int32_t value);
		void print(int16_t value);
		void print(int8_t value);

		char* read();

		uint16_t available();

};




#endif //SERIAL_H