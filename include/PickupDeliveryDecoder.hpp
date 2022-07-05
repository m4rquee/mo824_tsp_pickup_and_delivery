#ifndef PICKUPDELIVERYDECODER_DEFINE
#define PICKUPDELIVERYDECODER_DEFINE

#include "pickup_delivery_utils.hpp"
#include <algorithm>
#include <list>
#include <vector>

class PickupDeliveryDecoder {
    Pickup_Delivery_Instance &P;

public:
    explicit PickupDeliveryDecoder(Pickup_Delivery_Instance &P);
    ~PickupDeliveryDecoder();

    double decode(const std::vector<double> &chromosome) const;
    void decode(const std::vector<double> &chromosome, DNodeVector &Sol) const;
};

#endif// PICKUPDELIVERYDECODER_DEFINE
