//
// Created by tfuhrman on 09/05/17.
//

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "protocol.h"
#include <stdlib.h>
#include "canSender.h"
//get from old protocol file
// extern "C" {
    #include "compat.h"
    #include "robotstate.h"
    #include "control.h"
    #include "goals.h"
    #include "emergency.h"
//}

void autoSendStatus() {
    //todo auto_send ?
    CanSender::canSend(CURRENT_POS, 
        (int)current_pos.x, 
        (int)current_pos.y,
        (int)(current_pos.angle*FLOAT_PRECISION));
    
    CanSender::canSend(CURRENT_PWM, 
        (int)control.speeds.pwm_left, 
        (int)control.speeds.pwm_right);

    CanSender::canSend(CURRENT_SPD,
        (int)(control.speeds.linear_speed), 
        (int)wheels_spd.left, 
        (int)wheels_spd.right);

#if DEBUG_TARGET_SPEED
//    index += sprintf(message+index, ";%i;%i;%i;%i",
//			(int)wheels_spd.left,
//			(int)wheels_spd.right,
//			(int)(control.speeds.linear_speed - control.speeds.angular_speed),
//			(int)(control.speeds.linear_speed + control.speeds.angular_speed));
#endif
}


void ProtocolAutoSendStatus() {
    autoSendStatus();
// #if AUTO_STATUS_HZ
//     static int i=0;
//     if (++i % (HZ / AUTO_STATUS_HZ) == 0) {
//         autoSendStatus();
//     }
// #endif
}

uint8_t getLog10(const uint16_t number) {
    if(number>=10000) return 5;
    if(number>=1000) return 4;
    if(number>=100) return 3;
    if(number>=10) return 2;
    return 1;
}

unsigned char flagArduinoConnected = 0;

//order is order;id_servo;params
void parseAndExecuteOrder(uint8_t* message) {
    uint8_t mode = message[0];
    int order_id = 0;

    int8_t value;

    switch (mode) 
    {
        case HANDSHAKE:
        {
            // Ack that stm has started
            g_serial.print("HANDSHAKE");
            g_serial.print("\n");
            CanSender::canSend(WHOAMI,16);
            //g_serial.print("HANDSHAKE\n");
            // CanSender::canSend(SERIAL_INFO, "%d;", order_id);
            // CanSender::canSend(SERIAL_DEBUG, "Arduino %s has started (%d)", ARDUINO_ID, order_id);
            
            break;
        }

        // case WHOAMI:
        // {
        //     break;
        // }

        case SET_MODE:
        {
            break;
        }

        case SPEED:
        {
            g_serial.print("SPEED");
            g_serial.print("\n");
            int16_t l, a, t;
            //sscanf(receivedOrderPtr, "%i;%i;%i", &l, &a, &t);

            l =     message[1] << 8 |  message[2];

            a =     message[3] << 8 | message[4];

            t =     message[5] << 8 | message[6];

            g_serial.print(l);
            g_serial.print("\n");
            g_serial.print(a);
            g_serial.print("\n");
            g_serial.print(t);
            g_serial.print("\n");
            

            goal_data_t goal;
            goal.spd_data = {(float)t, l, a};
            FifoPushGoal(order_id, TYPE_SPD, goal);
            break;
        }

        case GET_CODER:
        {
            // CanSender::canSend(GET_CODER,get_left_encoder(), get_right_encoder());
            break;
        }

        case MANAGEMENT:
        {
            uint8_t order = message[1];
            switch(order)
            {
                case STOP:
                {
                    g_serial.print("STOP");
                    g_serial.print("\n");
                    CanSender::canSend(WHOAMI,CAN_ADDR);
                    flagArduinoConnected = 0;
                    break;
                }
                case START:
                {   
                    g_serial.print("START");
                    g_serial.print("\n");
                    CanSender::canSend(1,32);
                    flagArduinoConnected = 1;
                    break;
                }
                case PAUSE:
                {
                    g_serial.print("PAUSE");
                    g_serial.print("\n");
                    ControlSetStop(PAUSE_BIT);
                    break;
                }
                case RESUME:
                {
                    g_serial.print("RESUME");
                    g_serial.print("\n");
                    ControlUnsetStop(PAUSE_BIT);
                    break;
                }
                case RESET_ID:
                {
                    g_serial.print("RESET_ID");
                    g_serial.print("\n");
                    control.last_finished_id = 0;
                    break;
                }
                case SETEMERGENCYSTOP:
                {
                    g_serial.print("EMGSTOP");
                    g_serial.print("\n");
                    int enable = message[1];
                    EmergencySetStatus(enable);
                    break;
                }
                case NEXT_ORDER:
                {
                    g_serial.print("NEXT GOAL");
                    g_serial.print("\n");
                    FifoNextGoal();
                    ControlPrepareNewGoal();
                    break;
                }
                case RESET_ORDERS:
                {
                    g_serial.print("RESET GOALS");
                    g_serial.print("\n");
                    FifoClearGoals();
                    ControlPrepareNewGoal();
                    break;
                }
            }

            break;
        }

        case GOTOA:
        {
            g_serial.print("GOTOA");
            g_serial.print("\n");
            int16_t x, y, a_int, direction;
            float a;
            direction = 0;

            x =     message[1] << 8 | message[2];

            y =     message[3] << 8 | message[4];

            a_int =     message[5] << 8 | message[6];

            direction = message[7];

            g_serial.print(x);
            g_serial.print("\n");
            g_serial.print(y);
            g_serial.print("\n");
            g_serial.print(a_int);
            g_serial.print("\n");
            g_serial.print(direction);
            g_serial.print("\n");

            a = a_int / (float)FLOAT_PRECISION;
            goal_data_t goal;
            goal.pos_data = {x, y, direction};
            FifoPushGoal(order_id, TYPE_POS, goal);
            goal.ang_data = {a, 1};
            FifoPushGoal(order_id, TYPE_ANG, goal);
            break;
        }

        case GOTO:
        {
            g_serial.print("GOTO");
            g_serial.print("\n");
            int16_t x, y;
            int8_t direction;

            x =     message[1] << 8 | message[2];

            y =     message[3] << 8 | message[4];

            g_serial.print(x);
            g_serial.print("\n");
            g_serial.print(y);
            g_serial.print("\n");

            direction = message[5];
            
            goal_data_t goal;
            goal.pos_data = {x, y, direction};
            FifoPushGoal(order_id, TYPE_POS, goal);
            break;
        }

        case ROT:
        {
            g_serial.print("ROT");
            g_serial.print("\n");
            int16_t a_int;
            float a;

            a_int =     message[1] << 8 | message[2];

            g_serial.print(a_int);
            g_serial.print("\n");

            a = a_int / (float)FLOAT_PRECISION;
            goal_data_t goal;
            goal.ang_data = {a, 1};
            FifoPushGoal(order_id, TYPE_ANG, goal);
            break;
        }

        case ROTNOMODULO:
        {
            g_serial.print("ROT NO MODULO");
            g_serial.print("\n");
            long a_int;
            float a;
            
            a_int =     message[1] << 8 | message[2];
            g_serial.print(a_int);
            g_serial.print("\n");
            

            a = a_int / (float)FLOAT_PRECISION;
            goal_data_t goal;
            goal.ang_data = {a, 0};
            FifoPushGoal(order_id, TYPE_ANG, goal);
            break;
        }

        case PIDLEFT:
        case PIDRIGHT:
        case PIDALL:
        {
            long p_int, i_int, d_int;
            float p, i, d;
            
            p_int =     message[1] << 8 | message[2];

            i_int =     message[3] << 8 | message[4];

            d_int =     message[5] << 8 | message[6];

            g_serial.print(p_int);
            g_serial.print("\n");
            g_serial.print(i_int);
            g_serial.print("\n");
            g_serial.print(d_int);
            g_serial.print("\n");
            
            p = (float)p_int / (float)FLOAT_PRECISION;
            i = (float)i_int / (float)FLOAT_PRECISION;
            d = (float)d_int / (float)FLOAT_PRECISION;

            
            if (mode == PIDLEFT)
            {
                g_serial.print("PID LEFT");
                g_serial.print("\n");
                PIDSet(&PID_left, p, i, d, LEFT_BIAS);
            }
            else if (mode == PIDRIGHT)
            {
                g_serial.print("PID RIGHT");
                g_serial.print("\n");
                PIDSet(&PID_right, p, i, d, RIGHT_BIAS);
            }
            else 
            {
                g_serial.print("PID ALL");
                g_serial.print("\n");
                PIDSet(&PID_left, p, i, d, LEFT_BIAS);
                PIDSet(&PID_right, p, i, d, RIGHT_BIAS);
            }
            // CanSender::canSend(SERIAL_INFO, "%d;", order_id);
            break;
        }

        case PWM:
        {
            g_serial.print("PWM");
            g_serial.print("\n");
            int16_t l, r;
            uint16_t t;

            l =     message[1] << 8 | message[2];

            r =     message[3] << 8 | message[4];

            t =     message[5] << 8 | message[6];

            g_serial.print(l);
            g_serial.print("\n");
            g_serial.print(r);
            g_serial.print("\n");
            g_serial.print(t);
            g_serial.print("\n");

            goal_data_t goal;
            goal.pwm_data = {(float)t, l, r};
            FifoPushGoal(order_id, TYPE_PWM, goal);
            break;
        }

        case SET_POS:
        {
            g_serial.print("SET POS");
            g_serial.print("\n");
            int16_t x, y, a_int;
            float angle;

            
            x =     message[1] << 8 | message[2];

            y =     message[3] << 8 | message[4];

            a_int = message[5] << 8 | message[6];

            angle = a_int / (float)FLOAT_PRECISION;

            char buffer[200];
            int ret = snprintf(buffer, sizeof buffer, "%lf", angle);

            g_serial.print(x);
            g_serial.print("\n");
            g_serial.print(y);
            g_serial.print("\n");
            g_serial.print(a_int);
            g_serial.print("\n");

            



            RobotStateSetPos(x, y, angle);
            break;
        }

        case SET_PARAM:
        {
            g_serial.print("SET PARAM");
            g_serial.print("\n");
            int16_t r_int, s,a;
            float r;
            
            s =     message[1] << 8 | message[2];

            r_int = message[3] << 8 | message[4];

            a =     message[5] << 8 | message[6];

            g_serial.print(s);
            g_serial.print("\n");
            g_serial.print(r_int);
            g_serial.print("\n");
            g_serial.print(a);
            g_serial.print("\n");

            r = r_int / (float)FLOAT_PRECISION;
            control.max_spd = s;
            control.rot_spd_ratio = r;
            control.max_acc = a;
            break;
        }

        
        // case HALT:
        // {
        //     // Ack that arduino has stopped
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     CanSender::canSend(SERIAL_DEBUG, "Arduino %s has stopped (%d)", ARDUINO_ID, order_id);
        //     flagArduinoConnected = 0;
        //     break;
        // }
        // case PINGPING:
        //     //todo add LED on arduino
        //     digitalWrite(LED_DEBUG, HIGH);
        //     delay(1);
        //     digitalWrite(LED_DEBUG, LOW);
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     break;
        // case GET_CODER:
        //     CanSender::canSend(SERIAL_INFO, "%d;%d;%d;", order_id, left_ticks, right_ticks);
        //     break;
        // case GOTO:
        // {
        //     int x, y, direction;
        //     direction = 0;
        //     sscanf(receivedOrderPtr, "%i;%i;%i", &x, &y, &direction);
        //     goal_data_t goal;
        //     goal.pos_data = {x, y, direction};
        //     FifoPushGoal(order_id, TYPE_POS, goal);
        //     break;
        // }
        // case GOTOA:
        // {
        //     int x, y, a_int, direction;
        //     float a;
        //     direction = 0;
        //     sscanf(receivedOrderPtr, "%i;%i;%i;%i", &x, &y, &a_int, &direction);
        //     a = a_int / (float)FLOAT_PRECISION;
        //     goal_data_t goal;
        //     goal.pos_data = {x, y, direction};
        //     FifoPushGoal(order_id, TYPE_POS, goal);
        //     goal.ang_data = {a, 1};
        //     FifoPushGoal(order_id, TYPE_ANG, goal);
        //     break;
        // }
        // case ROT:
        // {
        //     int a_int;
        //     float a;
        //     sscanf(receivedOrderPtr, "%i", &a_int);
        //     a = a_int / (float)FLOAT_PRECISION;
        //     goal_data_t goal;
        //     goal.ang_data = {a, 1};
        //     FifoPushGoal(order_id, TYPE_ANG, goal);
        //     break;
        // }
        // case ROTNOMODULO:
        // {
        //     long a_int;
        //     float a;
        //     sscanf(receivedOrderPtr, "%li", &a_int);
        //     a = a_int / (float)FLOAT_PRECISION;
        //     goal_data_t goal;
        //     goal.ang_data = {a, 0};
        //     FifoPushGoal(order_id, TYPE_ANG, goal);
        //     break;
        // }
        // case PWM:
        // {
        //     int l, r, t;
        //     sscanf(receivedOrderPtr, "%i;%i;%i", &l, &r, &t);
        //     goal_data_t goal;
        //     goal.pwm_data = {(float)t, l, r};
        //     FifoPushGoal(order_id, TYPE_PWM, goal);
        //     break;
        // }
        // case SPD:
        // {
        //     int l, a, t;
        //     sscanf(receivedOrderPtr, "%i;%i;%i", &l, &a, &t);
        //     goal_data_t goal;
        //     goal.spd_data = {(float)t, l, a};
        //     FifoPushGoal(order_id, TYPE_SPD, goal);
        //     break;
        // }
        // case PIDALL:
        // case PIDRIGHT:
        // case PIDLEFT:
        // {
        //     long p_int, i_int, d_int;
        //     float p=0.5, i=0.5, d=0.5;
            
        //     sscanf(receivedOrderPtr, "%li;%li;%li", &p_int, &i_int, &d_int);
        //     // CanSender::canSend(SERIAL_INFO, "p: %l i: %l d: %l",p,i,d);
        //     // CanSender::canSend(SERIAL_INFO, "p: %(li) i: (%li) d: (%li)",p_int,i_int,d_int);
        //     p = (float)p_int / (float)FLOAT_PRECISION;
        //     i = (float)i_int / (float)FLOAT_PRECISION;
        //     d = (float)d_int / (float)FLOAT_PRECISION;
        //     //CanSender::canSend(SERIAL_INFO, "p: %f i: %f d: %f",p,i,d);
            
        //     if (orderChar == PIDLEFT)
        //         PIDSet(&PID_left, p, i, d, LEFT_BIAS);
        //     else if (orderChar == PIDRIGHT)
        //         PIDSet(&PID_right, p, i, d, RIGHT_BIAS);
        //     else {
        //         PIDSet(&PID_left, p, i, d, LEFT_BIAS);
        //         PIDSet(&PID_right, p, i, d, RIGHT_BIAS);
        //     }
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     Serial.println(p);
        //     Serial.println(i);
        //     Serial.println(d);
        //     break;
        // }
        // case KILLG:
        //     FifoNextGoal();
        //     ControlPrepareNewGoal();
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     break;
        // case CLEANG:
        //     FifoClearGoals();
        //     ControlPrepareNewGoal();
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     break;
        // case RESET_ID:
        //     control.last_finished_id = 0;
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     break;
        // case SET_POS:
        // {
        //     int x, y, a_int;
        //     float angle;
        //     sscanf(receivedOrderPtr, "%i;%i;%i;", &x, &y, &a_int);
        //     angle = a_int / (float)FLOAT_PRECISION;
        //     RobotStateSetPos(x, y, angle);
        //     CanSender::canSend(SERIAL_INFO, "%d;%i;", order_id, a_int);
        //     break;
        // }
        // case GET_POS:
        // {
        //     int x, y, a_int;
        //     float a;
        //     a = current_pos.angle;
        //     x = round(current_pos.x);
        //     y = round(current_pos.y);
        //     a_int = a * (float)FLOAT_PRECISION;
        //     CanSender::canSend(SERIAL_INFO, "%d;%d;%d;%d;", order_id, x, y, a_int);
        //     break;
        // }
        // case GET_SPD:
        // {
        //     int l, r;
        //     l = wheels_spd.left;
        //     r = wheels_spd.right;
        //     CanSender::canSend(SERIAL_INFO, "%d;%d;%d;", order_id, l, r);
        //     break;
        // }
        // case GET_TARGET_SPD:
        // {
        //     int left_spd, right_spd;
        //     left_spd = control.speeds.linear_speed - control.speeds.angular_speed;
        //     right_spd = control.speeds.linear_speed + control.speeds.angular_speed;
        //     CanSender::canSend(SERIAL_INFO, "%d;%d;%d;", order_id, left_spd, right_spd);
        //     break;
        // }
        // case GET_POS_ID:
        // {
        //     int x, y, a_int;
        //     float a;
        //     a = current_pos.angle;
        //     x = round(current_pos.x);
        //     y = round(current_pos.y);
        //     a_int = a * (float)FLOAT_PRECISION;
        //     CanSender::canSend(SERIAL_INFO, "%d;%d;%d;%d;", order_id, x, y, a_int, control.last_finished_id);
        //     break;
        // }
        // case SPDMAX:
        // {
        //     int r_int, s;
        //     float r;
        //     sscanf(receivedOrderPtr, "%i;%i", &s, &r_int);
        //     r = r_int / (float)FLOAT_PRECISION;
        //     control.max_spd = s;
        //     control.rot_spd_ratio = r;
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     break;
        // }
        // case ACCMAX:
        // {
        //     int a;
        //     sscanf(receivedOrderPtr, "%i", &a);
        //     control.max_acc = a;
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     break;
        // }
        // case GET_LAST_ID:
        //     CanSender::canSend(SERIAL_INFO, "%d", control.last_finished_id);
        //     break;
        // case PAUSE:
        //     ControlSetStop(PAUSE_BIT);
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     break;
        // case RESUME:
        //     ControlUnsetStop(PAUSE_BIT);
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     break;
        // case WHOAMI:
        //     CanSender::canSend(SERIAL_INFO, "%d;%s;", order_id, ARDUINO_ID);
        //     break;
        // case SETEMERGENCYSTOP:
        // {
        //     int enable;
        //     sscanf(receivedOrderPtr, "%i", &enable);
        //     EmergencySetStatus(enable);
        //     CanSender::canSend(SERIAL_INFO, "%d;", order_id);
        //     break;
        // }
        // default:
        //     CanSender::canSend(SERIAL_INFO, "Order %c is wrong !", orderChar);
    }
}


int8_t unsigned2signed(uint8_t number) 
{
    bool positive;
    int8_t result = 0;
    positive = number && MSB_ONLY_MASK == 0;
    // if (positive)
    // {
    //     result = MSB_POSITIVE;
    // }
    // else
    // {
    //     result = MSB_NEGATIVE;
    // }
    result = number;// && IGNORE_MSB_MASK;
    return result;
}