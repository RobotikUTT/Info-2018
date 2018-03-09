#ifndef __CAN_H__
#define __CAN_H__

/** Includes **/
/**************/
#include "stm32f3xx_hal.h"
#include "stm32f3xx_hal_can.h"
/** Defines **/
/*************/


/** Class Descritpion **/

class Can
{
	private:
		CAN_HandleTypeDef* m_can_interface_ptr;
		CanTxMsgTypeDef m_tx_msg;
		CanRxMsgTypeDef m_rx_msg;

	public:
		Can(CAN_HandleTypeDef* can, uint16_t id);
		~Can();

		uint8_t* read();
		void write(uint8_t* msg);
		uint8_t* get_rx_msg();
		uint8_t* get_tx_msg();

};


#endif