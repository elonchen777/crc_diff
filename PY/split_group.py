import pandas as pd


def prepare_early_late_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Create early and late groups based on tnm_stage and smoking_label"""
    def _label(row):
        crc = int(row.get('crc_label', 0))
        tnm_stage = row.get('tnm_stage', 0)
        smoking = int(row.get('smoking_label', 0))
        
        # Only consider CRC patients (crc_label == 1)
        if crc == 1:
            # Process TNM stage: early (stages 0,1) and late (stages 2,3,4)
            # Convert tnm_stage to integer for comparison
            try:
                tnm_int = int(tnm_stage)
            except (ValueError, TypeError):
                # Return None if cannot convert to integer
                return None
            
            if tnm_int in [0, 1]:
                stage_group = 'early'  # Early stage
            elif tnm_int in [2, 3, 4]:
                stage_group = 'late'   # Late stage
            else:
                # For other values, return None
                return None
            
            if smoking == 1:  # Smoking
                return f'CRC_{stage_group}_smoking'
            else:  # Non-smoking
                return f'CRC_{stage_group}_nonsmoking'
        else:
            # For non-CRC samples, return None (will be filtered later)
            return None
    
    df = df.copy()
    df['group'] = df.apply(_label, axis=1)
    # Filter out non-CRC samples and invalid TNM stages
    df = df[df['group'].notna()].copy()
    return df

def prepare_crc_smoke_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Create early and late groups based on tnm_stage and smoking_label.
    The age value is parsed and attached to each row (currently not used in the label).
    """
    def _label(row):
        crc = int(row.get('crc_label', 0))
        # tnm_stage = row.get('tnm_stage', 0)
        smoking = int(row.get('smoking_label', 0))
        # age = int(row.get('age', 0))
        
        # Only consider CRC patients (crc_label == 1)
        if crc == 1:
            if smoking == 1:  # Smoking
                return f'CRC_smoking'
            else:  # Non-smoking
                return f'CRC_nonsmoking'
        else:
            return f'CTRL'
    
    df = df.copy()
    df['group'] = df.apply(_label, axis=1)
    # Filter out non-CRC samples and invalid TNM stages
    df = df[df['group'].notna()].copy()
    return df

def prepare_crc_diff_groups(df: pd.DataFrame) -> pd.DataFrame:
    def _label(row):
        crc = int(row.get('crc_label', 0))
        diff = int(row.get('differentiation', 0))
        
        # Only consider CRC patients (crc_label == 1)
        if crc == 1:
            if diff == 1:
                return f'CRC_poordiff'
            else:
                return f'CRC_welldiff'
        else:
            return f'CTRL'
    
    df = df.copy()
    df['group'] = df.apply(_label, axis=1)
    df = df[df['group'].notna()].copy()
    return df