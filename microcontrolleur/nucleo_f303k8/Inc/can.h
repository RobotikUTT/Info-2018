#ifndef __CAN_H__
#define __CAN_H__

/** Includes **/
/**************/
#include "stm32f3xx_hal.h"
#include "stm32f3xx_hal_can.h"
/** Defines **/
/*************/


/** Class Descritpion **/
extern uint16_t g_nb_msg_received;
class Can
{
	private:
		CAN_HandleTypeDef* m_can_interface_ptr;
		CanTxMsgTypeDef m_tx_msg;
		CanRxMsgTypeDef m_rx_msg;
		uint8_t prev_msg[8]= {0,0,0,0,0,0,0,0};

		

		bool check_eq_msgs(uint8_t* msg1, uint8_t* msg2,uint8_t size);
		void copy_msg(uint8_t* dest,uint8_t* src, uint8_t size);

	public:
		Can(CAN_HandleTypeDef* can, uint16_t id);
		~Can();

		uint8_t* read();
		HAL_StatusTypeDef write(uint8_t* msg);
		// uint8_t* get_rx_msg();
		uint16_t available();

};


#endif