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
	// m_rx_msg.IDE = CAN_ID_STD;
	// m_rx_msg.RTR = CAN_RTR_DATA;
	// m_rx_msg.FMI = 2;
	m_can_interface_ptr->pTxMsg = &m_tx_msg;
	m_can_interface_ptr->pRxMsg = &m_rx_msg;

}

Can::~Can()
{
	delete m_can_interface_ptr;
}
/** Public Methods **/
/********************/

void Can::write(uint8_t* msg)
{
 

  for( int i = 0; i < 8; i++)
  {
  	m_tx_msg.Data[i] = msg[i];
  } 
  
  // hcan.pTxMsg = &m_tx_msg;

  HAL_StatusTypeDef can_status;
  // can_status = HAL_CAN_Transmit(m_can_interface_ptr,0x0FFF);
  
  switch(can_status)
  {


    case HAL_OK:
	    // g_serial.print("CAN SENT OK\n");
	    // HAL_UART_Transmit(&huart2, "CAN SEND OK", 12, 0xFFFF);
    break;

    case HAL_ERROR:
    // strcpy(msg, "CAN SENT ERROR\nERROR CODE:" );
    	// HAL_UART_Transmit(&huart2, "CAN SEND ER", 12, 0xFFFF);
    //send_PC_UART_DATA((long)HAL_CAN_GetError(&hcan));
    // send_uint32(HAL_CAN_GetError(&hcan));
    // strcpy(msg, "\n" );
    // HAL_UART_Transmit(&huart2, (uint8_t*)msg, strlen(msg), 0xFFFF);
    	break;

    case HAL_BUSY:
    	// HAL_UART_Transmit(&huart2, "CAN SEND BS", 12, 0xFFFF);
    // strcpy(msg, "CAN BUSY\n");
    // HAL_UART_Transmit(&huart2, (uint8_t*)msg, strlen(msg), 0xFFFF);
    	break;

    case HAL_TIMEOUT:
    	// HAL_UART_Transmit(&huart2, "CAN SEND TO", 12, 0xFFFF);
    // strcpy(msg, "CAN TIMEOUT\n");
    // HAL_UART_Transmit(&huart2, (uint8_t*)msg, strlen(msg), 0xFFFF);
    	break;

  }
}

uint8_t* Can::read()
{

  // HAL_StatusTypeDef can_status;
  // can_status = HAL_CAN_Receive(m_can_interface_ptr, 50,0x0FFF);
  // switch(can_status)
  // {
  //   case HAL_OK:
  //   strcpy(msg, "CAN RECEIVED OK\n");
  //   m_rx_msg = *hcan.pRxMsg;
  //   HAL_UART_Transmit(&huart2, (uint8_t*)msg, strlen(msg), 0xFFFF);
  //   break;

  //   case HAL_ERROR:
  //   strcpy(msg, "CAN RECEIVED ERROR\nERROR CODE:" );
  //   HAL_UART_Transmit(&huart2, (uint8_t*)msg, strlen(msg), 0xFFFF);
  //   //send_PC_UART_DATA((long)HAL_CAN_GetError(&hcan));
  //   send_uint32(HAL_CAN_GetError(&hcan));
  //   strcpy(msg, "\n" );
  //   HAL_UART_Transmit(&huart2, (uint8_t*)msg, strlen(msg), 0xFFFF);
  //   break;

  //   case HAL_BUSY:
  //   strcpy(msg, "CAN BUSY\n");
  //   HAL_UART_Transmit(&huart2, (uint8_t*)msg, strlen(msg), 0xFFFF);
  //   break;

  //   case HAL_TIMEOUT:
  //   strcpy(msg, "CAN TIMEOUT\n");
  //   HAL_UART_Transmit(&huart2, (uint8_t*)msg, strlen(msg), 0xFFFF);
  //   break;

  //   default:
  //     strcpy(msg, "CAN RECEPTION DEFAULT\n");
  //     HAL_UART_Transmit(&huart2, (uint8_t*)msg, strlen(msg), 0xFFFF);
  //     break;

  // }
  return 0;// (m_rx_msg->Data);
  
}

uint8_t* Can::get_rx_msg()
{
	return m_rx_msg.Data;
}

uint8_t* Can::get_tx_msg()
{
	return m_tx_msg.Data;
}

/** Private Methods **/
/*********************/