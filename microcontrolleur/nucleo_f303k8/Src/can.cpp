/** Includes **/
/**************/
#include "can.h"

/** Constructor **/
/*****************/

Can::Can(CAN_HandleTypeDef* can, uint16_t id)
{
	m_can_interface_ptr = can;

	// tx message init
	m_tx_msg.StdId = id;
  m_tx_msg.ExtId = id;
	m_tx_msg.IDE = CAN_ID_STD;
	m_tx_msg.RTR = CAN_RTR_DATA;
	m_tx_msg.DLC = 8;
	
	// rx message init  	
	m_rx_msg.DLC = 8;
	m_rx_msg.IDE = CAN_ID_STD;
  m_rx_msg.RTR = CAN_RTR_DATA;
  m_rx_msg.FMI = 0;
  m_rx_msg.FIFONumber = CAN_FIFO0;

	m_can_interface_ptr->pTxMsg = &m_tx_msg;
	m_can_interface_ptr->pRxMsg = &m_rx_msg;

}

Can::~Can()
{
	delete m_can_interface_ptr;
}
/** Public Methods **/
/********************/

HAL_StatusTypeDef Can::write(uint8_t* msg)
{
 
  for( int i = 0; i < 8; i++)
  {
  	m_tx_msg.Data[i] = msg[i];
  } 
  
  return HAL_CAN_Transmit_IT(m_can_interface_ptr);
  
}

HAL_StatusTypeDef Can::read()
{

  
  HAL_StatusTypeDef _can_status;
  _can_status = HAL_CAN_Receive_IT(m_can_interface_ptr, CAN_FIFO0);

  return _can_status;
  
}

uint8_t* Can::get_rx_msg()
{
	return m_rx_msg.Data;
}

/** Private Methods **/
/*********************/