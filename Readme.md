#Parallel Monocular SLAM

## Optimisations

### Bundle Adjustment 
``optimiseInnerWindow``:
- Uses keyframes and frames
- Optimises both points reconstructed in each keyframe and the 3d features. These are apparently different. It looks like features are turned into points when they are matched between keyframes. Features are added when there are matched measurements between the keyframe and its small baseline frame.

